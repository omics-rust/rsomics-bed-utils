//! Group tab-delimited rows by key column(s) and aggregate a value column.
//!
//! Matches `bedtools groupby`: consecutive rows with identical key-column values
//! form a group; when the key changes the group is flushed. Rows must therefore
//! be pre-sorted on the group key (same contract as bedtools).

#![allow(clippy::cast_precision_loss)]

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write, stdin};

use rsomics_common::{Result, RsomicsError};

/// Parse a comma-or-range list of 1-based column indices like "1,2" or "1-3".
fn parse_col_list(s: &str) -> Result<Vec<usize>> {
    let mut cols = Vec::new();
    for part in s.split(',') {
        if let Some((lo, hi)) = part.split_once('-') {
            let lo: usize = lo
                .trim()
                .parse()
                .map_err(|_| RsomicsError::InvalidInput(format!("bad column range: {part:?}")))?;
            let hi: usize = hi
                .trim()
                .parse()
                .map_err(|_| RsomicsError::InvalidInput(format!("bad column range: {part:?}")))?;
            if lo == 0 || hi == 0 || lo > hi {
                return Err(RsomicsError::InvalidInput(format!(
                    "invalid column range: {part:?}"
                )));
            }
            cols.extend(lo..=hi);
        } else {
            let n: usize = part
                .trim()
                .parse()
                .map_err(|_| RsomicsError::InvalidInput(format!("bad column index: {part:?}")))?;
            if n == 0 {
                return Err(RsomicsError::InvalidInput(
                    "column indices are 1-based".to_string(),
                ));
            }
            cols.push(n);
        }
    }
    Ok(cols)
}

/// Format a float matching bedtools' `%.10g` output: integer-valued results
/// print without a decimal point; non-integer values use 10 significant figures
/// with trailing zeros stripped.
fn fmt_float(v: f64) -> String {
    if v.fract() == 0.0 && v.abs() < 1e15 {
        return format!("{}", v as i64);
    }
    // Compute the number of decimal places needed for 10 significant figures.
    let magnitude = if v == 0.0 {
        0i32
    } else {
        v.abs().log10().floor() as i32
    };
    let decimals = (9 - magnitude).max(0) as usize;
    format!("{:.prec$}", v, prec = decimals)
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

fn aggregate(values: &[f64], op: &str) -> String {
    match op {
        "sum" => fmt_float(values.iter().sum()),
        "mean" => {
            if values.is_empty() {
                return ".".to_string();
            }
            fmt_float(values.iter().sum::<f64>() / values.len() as f64)
        }
        "min" => {
            if values.is_empty() {
                return ".".to_string();
            }
            fmt_float(values.iter().copied().fold(f64::INFINITY, f64::min))
        }
        "max" => {
            if values.is_empty() {
                return ".".to_string();
            }
            fmt_float(values.iter().copied().fold(f64::NEG_INFINITY, f64::max))
        }
        "count" => format!("{}", values.len()),
        "stdev" => {
            if values.len() < 2 {
                return "0".to_string();
            }
            let mean = values.iter().sum::<f64>() / values.len() as f64;
            let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64;
            fmt_float(var.sqrt())
        }
        _ => ".".to_string(),
    }
}

fn aggregate_str(raw: &[String], op: &str) -> String {
    match op {
        "collapse" => raw.join(","),
        "distinct" => {
            let mut seen: HashSet<&str> = HashSet::new();
            let mut result = Vec::new();
            for s in raw {
                if seen.insert(s.as_str()) {
                    result.push(s.as_str());
                }
            }
            result.join(",")
        }
        "count_distinct" => {
            let set: HashSet<&str> = raw.iter().map(|s| s.as_str()).collect();
            format!("{}", set.len())
        }
        "first" => raw.first().map(|s| s.as_str()).unwrap_or(".").to_string(),
        "last" => raw.last().map(|s| s.as_str()).unwrap_or(".").to_string(),
        _ => ".".to_string(),
    }
}

fn is_numeric_op(op: &str) -> bool {
    matches!(op, "sum" | "mean" | "min" | "max" | "count" | "stdev")
}

struct Group {
    key_fields: Vec<String>,
    /// raw string values for each (col, op) pair
    raw_values: Vec<Vec<String>>,
    /// parsed float values (for numeric ops)
    num_values: Vec<Vec<f64>>,
}

impl Group {
    fn new(key_fields: Vec<String>, n_ops: usize) -> Self {
        Self {
            key_fields,
            raw_values: vec![Vec::new(); n_ops],
            num_values: vec![Vec::new(); n_ops],
        }
    }

    fn flush(&self, col_indices: &[usize], ops: &[String], out: &mut impl Write) -> Result<()> {
        let keys = self.key_fields.join("\t");
        let mut agg_parts = Vec::with_capacity(ops.len());
        for (i, op) in ops.iter().enumerate() {
            let part = if is_numeric_op(op) {
                aggregate(&self.num_values[i], op)
            } else {
                aggregate_str(&self.raw_values[i], op)
            };
            agg_parts.push(part);
        }
        // Suppress unused-variable warning; col_indices validated at parse time.
        let _ = col_indices;
        writeln!(out, "{keys}\t{}", agg_parts.join("\t")).map_err(RsomicsError::Io)
    }
}

fn groupby_inner(
    lines: impl Iterator<Item = std::io::Result<String>>,
    group_cols: &[usize],
    col_indices: &[usize],
    ops: &[String],
    out: &mut dyn Write,
) -> Result<()> {
    let mut writer = BufWriter::with_capacity(64 * 1024, out);
    let mut current: Option<Group> = None;

    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();

        // Build this row's group key.
        let key: Vec<String> = group_cols
            .iter()
            .map(|&c| {
                fields
                    .get(c.saturating_sub(1))
                    .copied()
                    .unwrap_or("")
                    .to_string()
            })
            .collect();

        // If key changed, flush previous group.
        if let Some(ref g) = current
            && g.key_fields != key
        {
            g.flush(col_indices, ops, &mut writer)?;
            current = None;
        }

        let group = current.get_or_insert_with(|| Group::new(key.clone(), ops.len()));

        for (i, &col) in col_indices.iter().enumerate() {
            let raw = fields
                .get(col.saturating_sub(1))
                .copied()
                .unwrap_or("")
                .to_string();
            group.raw_values[i].push(raw.clone());
            if is_numeric_op(&ops[i])
                && let Ok(v) = raw.parse::<f64>()
            {
                group.num_values[i].push(v);
            }
        }
    }

    if let Some(g) = current {
        g.flush(col_indices, ops, &mut writer)?;
    }

    writer.flush().map_err(RsomicsError::Io)
}

pub fn groupby(
    input: &str,
    output: &mut dyn Write,
    group_spec: &str,
    col_spec: &str,
    op_spec: &str,
) -> Result<()> {
    let group_cols = parse_col_list(group_spec)?;
    let col_indices = parse_col_list(col_spec)?;
    let ops: Vec<String> = op_spec.split(',').map(|s| s.trim().to_string()).collect();

    if col_indices.len() != ops.len() {
        return Err(RsomicsError::InvalidInput(
            "-c and -o must have the same number of comma-separated entries".to_string(),
        ));
    }

    if input == "-" {
        groupby_inner(
            stdin().lock().lines(),
            &group_cols,
            &col_indices,
            &ops,
            output,
        )
    } else {
        let file =
            File::open(input).map_err(|e| RsomicsError::InvalidInput(format!("{input}: {e}")))?;
        groupby_inner(
            BufReader::new(file).lines(),
            &group_cols,
            &col_indices,
            &ops,
            output,
        )
    }
}
