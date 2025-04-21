use std::path::PathBuf;
pub mod vcf;
use tokio;

#[tokio::main]
async fn main() {
    let path = PathBuf::from(
        r"/Users/jss009/code/vrs_anvil_toolkit/u08_release_data/gregor_consortium_u06_sorted_chr1_V2_VT_VEP_VRS.vcf.gz",
    );
    let _ = vcf::load_vcf(path, None, PathBuf::from("output.txt")).await;
    ()
}
