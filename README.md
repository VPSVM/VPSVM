# Introduction
This codebase is part of the VPSVM paper.

**What this codebase includes**: example and benchmark implementations in C++17 for some of the schemes in the VPSVM paper.

**What this codebase is not**: it is not for production use; it is not extensively tested.

# Platforms
In our paper, we limit the frequency of the CPU to 800MHZ, in order to simulate the resource-restricted users. If you want to get the same test results as in the paper, you need to limit the CPU frequency to 800MHZ when executing the algorithm on the user side.

# IDE  
In the coding, we chose CLion as IDE. Therefore, if you use CLion as IDE it will help you deploy this project more conveniently.

# Dependencies  
## Install m4
    wget http://mirrors.kernel.org/gnu/m4/m4-1.4.13.tar.gz
    tar -xzvf m4-1.4.13.tar.gz
    cd m4-1.4.13
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
## Install autoconf
    wget http://mirrors.kernel.org/gnu/autoconf/autoconf-2.65.tar.gz
    tar -xzvf autoconf-2.65.tar.gz
    cd autoconf-2.65
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
## Install automake
    wget http://mirrors.kernel.org/gnu/automake/automake-1.11.tar.gz
    tar xzvf automake-1.11.tar.gz
    cd automake-1.11
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
## Install libtool
    wget http://mirrors.kernel.org/gnu/libtool/libtool-2.2.6b.tar.gz
    tar xzvf libtool-2.2.6b.tar.gz
    cd libtool-2.2.6b
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
## Install GMP
    tar -jxvf  gmp-6.2.1.tar.bz2
    cd gmp-6.2.1
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
    make check
## Install CMake
    tar -xzvf cmake-3.10.2.tar.gz
    cd cmake-3.10.2
    ./configure --prefix=/usr/local
    sudo make
    sudo make install
## Install NTL
    tar -xzvf ntl-11.5.1.tar.gz
    cd ntl-11.5.1/src
    ./configure 
    make
    sudo make install



