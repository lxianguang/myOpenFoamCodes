## Table of contents
  * [The weight pressure source theory](#WPS-Theory)
  * [Libtorch for OpenFOAM v2406](#Libtorch)


## WPS Theory
 to do !

## Libtorch
LibTorch can be used in the environment of OpenFOAM. Instead of using PyTorch, using libtorch seems to be straightforward in the environment of OpenFOAM since both of them are written in C++. In this tutorial, the programming of LibTorch and OpenFOAM is demonstrated step by step. It was tested by OpenFOAM v2406.

1. Download libtorch: https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-2.9.1%2Bcpu.zip
(unzip it, in my case, I put it under: my/path/OpenFOAM/libtorch)：
    ```
    unzip libtorch-shared-with-deps-2.9.1+cpu.zip
    ```

2. download this test case, myLibtorchFoam. (In this solver, you will see how to set Make/options in order to let OpenFOAM find libtorch)
    ```
    cd your/path/myLibtorchFoam
    wmake
    ```

3. wmake torchFoam.
    Open .bashrc file, add the following:
    ```
    export LD_LIBRARY_PATH=~/OpenFOAM/libtorch/lib:$LD_LIBRARY_PATH
    source ~/.bashrc
    ```

4. run myLibtorchFoam, output：
    ```
    0.2870  0.5473  0.5788
    0.5582  0.2020  0.7702
    [ CPUFloatType{2,3} ] 
    ```
