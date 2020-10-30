# EuroPLEx School 2020

## High Performance Computing and Aspects of Computing in Lattice Gauge Theories

### Preparation for practicals

#### Prerequisites 
You will need [`make`](https://www.gnu.org/software/make/), a C compiler, a Python (version 3) interpreter, and the python packages [`matplotlib`](https://matplotlib.org) and [`numpy`](https://matplotlib.org). The C compiler should be able to parse [OpenMP](https://www.openmp.org) pragmas; for the GNU C compiler, this means it should understand the `-fopenmp` option. Also have [Git](https://git-scm.com) installed, to be able to clone this repository, both to check your installation and to download the exercises during the lectures.

##### Linux
If on Linux, check your distribution's documentation and package manager for installing the above packages. For example, on Ubuntu:
```bash
[user@localhost ~]$ sudo apt install build-essential python3-numpy python3-matplotlib git-all
```
should be sufficient.

##### Windows
On Windows, you can use the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) (WSL). Follow [these instructions](https://docs.microsoft.com/en-us/windows/wsl/install-win10) and on Step 6. choose "Ubuntu 20.04 LTS". Follow the link to [Create a user account and password for your new Linux distribution](https://docs.microsoft.com/en-us/windows/wsl/user-support). Once you have run `sudo apt update` and  `sudo apt upgrade` then install the packages listed above as for the case of Linux.

##### MacOS
On MacOS follow the instructions for [installing Homebrew](https://docs.brew.sh/Installation), which includes installing the Command Line Tools for Xcode. The default LLVM clang compiler may not include OpenMP support. In this case, use Homebrew to install `libomp`. You can then use either homebrew or [pip](https://pip.pypa.io/en/stable/installing/) to install the required Python packages.

#### Test your installation
Clone this repo and change into it:
```bash
[user@localhost ~]$ git clone https://github.com/g-koutsou/EuroPLEx2020
Cloning into 'EuroPLEx2020'...
remote: Enumerating objects: 6, done.
remote: Counting objects: 100% (6/6), done.
remote: Compressing objects: 100% (6/6), done.
remote: Total 6 (delta 0), reused 6 (delta 0), pack-reused 0
Unpacking objects: 100% (6/6), 2.19 KiB | 2.19 MiB/s, done.
[user@localhost ~]$ cd EuroPLEx2020/
```
To test your C installation, build the example:
```bash
[user@localhost EuroPLEx2020]$ make
gcc -Xpreprocessor -fopenmp -c runme.c
gcc -o runme runme.o -fopenmp
[user@localhost EuroPLEx2020]$ 
```
If this fails, you may need to adjust compiler flags in the Makefile according to the compiler you are using.

Run using two OpenMP threads. You should see output as follows (the order of the output may differ):
```bash
[user@localhost EuroPLEx2020]$ OMP_NUM_THREADS=2 ./runme
n_threads = 2 thread_id = 0
n_threads = 2 thread_id = 1
```

To test your Python installation, run the python script:
```bash
[user@localhost EuroPLEx2020]$ python3 test-plot.py
```
This should generate a PDF plot `test-plot.pdf` in the current directory:
```bash
[user@localhost EuroPLEx2020]$ ls -1
Makefile
README.md
runme
runme.c
runme.o
test-plot.pdf
test-plot.py
```
The plot should look like this (typeface and color may differ):
![test-plot.pdf](https://i.imgur.com/glh90x4.png)
