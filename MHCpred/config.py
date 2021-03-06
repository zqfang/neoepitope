
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("--data_path", type=str, required=True)
parser.add_argument("--mhc", type=str, help="mapping of mhc allele name to pseudo_seq")
parser.add_argument("--log_dir", type=str, default="checkpoints")
parser.add_argument("--batch_size", type=int, default= 100)
args = parser.parse_args()


## Training settings
DATA_BUNDLE = args.data_path
# DATA_BUNDLE = "/data/bases/fangzq/ImmunoRep/IEDB/"
# mhc_allel_filename = "IEDB/data/allelenames"
# pesudo_filename = "IEDB/data/MHC_pseudo.dat"
# peptide_ba_el_dir = "IEDB/data/threshold/"
mhc2psedo_filename = args.mhc #"/data/bases/fangzq/ImmunoRep/IEDB/mhc2pseudo.csv"
eval_filename = args.data_path #"/data/bases/fangzq/ImmunoRep/IEDB/BD2013.bi.humanHLA.txt"
## Hyperparameters
class_num = 1
input_size = 1900
batch_size = args.batch_size
num_workers = 5
learning_rate = 0.001
num_epochs = 100
steps_per_validation  = 500


