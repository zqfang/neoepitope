

# PointNovo

This is a fork of [PointNovo](https://github.com/volpato30/PointNovo).   

**NOTE**: This repo will be refactor for neoantigen discovery.

Please refer to PointNovo to get details.



## Dependency
- python >= 3.7
- pytorch >= 1.0
- biopython
- pyteomics
- cython

For database search you also need to install [percolator](http://percolator.ms/).


## knapsack files
Like DeepNovo, in PointNovo we also use the knapsack algorithm to further limit the search space. This means when performing de novo sequencing,
the program needs to either read or create a knapsack matrix based on the selected PTMs (one time computation). Pre-built knapsack matrix files could be found [here](https://1drv.ms/u/s!AvnYi33QHIzqwyaJdF89AneoTVUY?e=BJCHqZ):

You can use symbolic links to choose which knapsack file to use. i.e.

~~~
ln -s fix_C_var_NMQ_knapsack.npy knapsack.npy
~~~

## usage
### first build cython modules

~~~
python deepnovo_cython_setup.py build_ext --inplace
~~~

### train mode:

~~~
python main.py --train
~~~

On a RTX 2080 Ti GPU it takes around 0.4 second to train a batch of 16 annotated spectra

### denovo mode:

~~~
python main.py --search_denovo
~~~

### evaluate denovo result:

~~~
python main.py --test
~~~

This script is borrowed from the original DeepNovo implementation. It will generate the metrics defined by the paper.

### database search mode:

~~~
python main.py --search_db
~~~




