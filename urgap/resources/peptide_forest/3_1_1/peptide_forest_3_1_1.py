"""PeptideForest resource."""

import argparse
import multiprocessing as mp

import peptide-forest 

if __name__ == "__main__":
    mp.freeze_support()
    mp.set_start_method("fork", force=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="config", help="config path json")
    parser.add_argument("-o", dest="output", help="output file")
    args = parser.parse_args()

    pf = peptide-forest.PeptideForest(
        config_path=args.config,
        output=args.output,
    )
    pf.prep_ursgal_csvs()
    pf.calc_features()
    pf.fit()
    pf.get_results()
    pf.write_output()
