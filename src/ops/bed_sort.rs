use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::bed;

pub fn sort_bed(input: &Path, output: &mut dyn Write) -> Result<()> {
    bed::sort_bed3(File::open(input).map_err(RsomicsError::Io)?, output)
}

pub fn sort_bed_stdin(output: &mut dyn Write) -> Result<()> {
    bed::sort_bed3(io::stdin().lock(), output)
}
