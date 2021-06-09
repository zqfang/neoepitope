import argparse
import logging
import logging.config
from train_func import train

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spectrum", required=True, help="mgf file")
    parser.add_argument("--train_feature", required=True, help="features with mass corrected csv file.")
    parser.add_argument("--valid_feature", required=True, help="features with mass corrected csv file.")
    parser.add_argument("--location_dict", required=True, help="feature's specum locaion of the mf file." )
    parser.add_argument("--knapsack", required=False, default="knapsack.npy", help="use the knapsack algorithm to limit the search space." )

    args = parser.parse_args()
    ## start training
    train(args)

if __name__ == '__main__':
    log_file_name = 'PointNovo.train.log'
    d = {
        'version': 1,
        'disable_existing_loggers': False,  # this fixes the problem
        'formatters': {
            'standard': {
                'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
            },
        },
        'handlers': {
            'console': {
                'level': 'INFO',
                'class': 'logging.StreamHandler',
                'formatter': 'standard',
            },
            'file': {
                'level': 'DEBUG',
                'class': 'logging.FileHandler',
                'filename': log_file_name,
                'mode': 'w',
                'formatter': 'standard',
            }
        },
        'root': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
        }
    }
    logging.config.dictConfig(d)
    main()