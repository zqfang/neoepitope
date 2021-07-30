
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("--data_bundle", type=str, required=True)
parser.add_argument("--log_dir", type=str, default="checkpoints")
parser.add_argument("--train", dest="train", action="store_true")
parser.add_argument("--batch_size", type=int, default= 1000)
parser.add_argument("--eval", dest="eval", action="store_true")
args = parser.parse_args()


## Training settings
DATA_BUNDLE = args.data_bundle
# DATA_BUNDLE = "/data/bases/fangzq/ImmunoRep/IEDB/databundle.pkl"
# mhc_allel_filename = "IEDB/data/allelenames"
# pesudo_filename = "IEDB/data/MHC_pseudo.dat"
# peptide_ba_el_dir = "IEDB/data/threshold/"

## Hyperparameters
class_num = 1
input_size = 1900
batch_size = args.batch_size
num_workers = 5
learning_rate = 0.001
num_epochs = 100


