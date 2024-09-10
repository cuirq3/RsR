<div align="center">

# RsR: Rotation-System Reconstruction

This is the official implementation of

## Surface Reconstruction Using Rotation Systems

<font size="4">
<a href="https://cuirq3.github.io/" style="font-size:100%;">Ruiqi Cui<sup>1</sup></a>&emsp;
<a href="https://orbit.dtu.dk/en/persons/emil-toftegaard-g%C3%A6de" style="font-size:100%;">Emil Toftegaard Gæde<sup>1</sup> </a>&emsp;
<a href="http://www2.compute.dtu.dk/~erot/" style="font-size:100%;">Eva Rotenberg<sup>1</sup> </a>&emsp;
<a href="https://www.graphics.rwth-aachen.de/person/3/" style="font-size:100%;">Leif Kobbelt<sup>2</sup> </a>&emsp;
<a href="https://people.compute.dtu.dk/janba/" style="font-size:100%;">J. Andreas Bærentzen<sup>1</sup> </a>&emsp;
</font>
<br>

<font size="4">
<sup>1</sup> Department of Applied Mathematics and Computer Science, Technical University of Denmark, Denmark

<sup>2</sup> Visual Computing Institute, RWTH Aachen University, Germany
</font>

| <a href="">Webpage(TODO)</a> | <a href="https://arxiv.org/abs/2402.01893">arXiv</a> | <a href="">Presentation video(TODO, link)</a> |

<img src="./pics/teaser.png" alt="teaser.png"/><img src="./pics/Scene_Stanford.png" alt="stanford.png"/> <br>
<!-- <b>Our method extracts meshes from 3D Gaussian Splatting reconstructions and builds hybrid representations <br>that enable easy composition and animation in Gaussian Splatting scenes by manipulating the mesh.</b> -->
</div>

# Abstract

_We propose a combinatorial method for reconstruction of surfaces from points. Our method
constructs a spanning tree and a rotation system. Since the tree is trivially
a planar graph, its rotation system determines a genus zero surface with a
single face which we proceed to incrementally refine by inserting edges to
split faces. In order to raise the genus, special handles are added in a later
stage by inserting edges between different faces and thus merging them. It turns out that our approach has two specific benefits over the other methods. First, the output mesh preserves the
most information from the input point cloud. Second, our method provides
control over the topology of the reconstructed surface._

# Updates and To-do list

<details>
<summary><span style="font-weight: bold;">Updates</span></summary>
<ul>
  <li><b>[09/09/2024]</b> Code release.</li>
</ul>
</details><br>

<details>
<summary><span style="font-weight: bold;">To-do list</span></summary>
<ul>
  <li><b>Improvement:</b> Remove the reliance on 3rd party libraries.</li>
</ul>
</details>
<br>

# Getting Started
We tested our project on [Windows](#Windows), [Linux](#linux), and [Mac](#mac) OS. You can find the corresponding installation guide below.

<details>

<summary><strong style="font-size: 1.5em;">&nbsp;Windows</strong></summary>

We recommend using <a href="https://github.com/microsoft/vcpkg?tab=readme-ov-file">VCPKG</a> + <a href="https://cmake.org/">CMake</a> for installation. A good tutorial can be found <a href="https://learn.microsoft.com/vcpkg/get_started/get-started">here</a>. We specify the software and library versions used in our tests, but users are not limited to these versions.

### 0. Prerequisites

All the libs are installed via VCPKG.

- CMake 3.29.0-rc2
- Visual Studio 2019
- VCPKG 
- Libs
  - [CGAL](https://www.cgal.org/download/windows.html) 5.6
  - [Boost](https://www.boost.org/doc/libs/1_86_0/more/getting_started/windows.html)
  - [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page) 3.4.0

### 1. Installing

- Clone the repo.
```
git clone https://github.com/cuirq3/RsR.git
cd RsR
```

- Configure and Generate - replace `path/to/vcpkg` with your local setting, e.g. `C:/vcpkg`
```
mkdir build
cd build
cmake -DCMAKE_TOOLCHAIN_FILE=path/to/vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Release ..
```
- Build - Manually do it in the IDE or run the following command
```
cmake --build . --config Release
```
- Done! You can jump to [A simple example](#a-simple-example) to test if it is successfully installed.
</details>


<details>

<summary><strong style="font-size: 1.5em;">&nbsp;Linux (Ubuntu)</strong></summary>

### 0. Prerequisites

Install required software and packages:

```
sudo apt install cmake
sudo apt install libcgal-dev
sudo apt install libboost-all-dev
sudo apt install libeigen3-dev
```

### 1. Installing

- Clone the repo.
```
git clone https://github.com/cuirq3/RsR.git
cd RsR
```

- Compile
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 12
```

- Done! You can go to [A simple example](#a-simple-example) to test if it is successfully installed.

</details>

<details>

<summary><strong style="font-size: 1.5em;">&nbsp;Mac</strong></summary>

We tested with <a href="https://www.macports.org">MacPorts</a> + <a href="https://cmake.org/">CMake</a> for installation. We specify the software and library versions used in our tests, but users are not limited to these versions.

### 0. Prerequisites

Install required software and packages:

```
sudo port install cmake
sudo port install cgal5
sudo port install boost
sudo port install eigen3
```

### 1. Installing

- Clone the repo.
```
git clone https://github.com/cuirq3/RsR.git
cd RsR
```

- Compile
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 12
```

- Done! You can go to [A simple example](#a-simple-example) to test if it is successfully installed.

</details>
<br>

# A simple example
After installation, a great test would be to reconstruct the well-known [Stanford Bunny](https://graphics.stanford.edu/data/3Dscanrep/). The point cloud is already included in the example folder. You can directly start reconstructing by running the following command.

Note: the commands should follow right after the installation. (Make sure you are at the correct working directory)
```
// For Windows
./Release/RsR.exe ../configs/example_config.txt

// For Linux and Mac
./RSR ../configs/example_config.txt
```

A reconstructed bunny should be found in the **output** folder.

# Config file
You can modify the config file to adapt the software to different scenarios.

# Results reproduction

# Acknowledgments

# BibTeX

```
@article{cui2024surface,
  title={Surface Reconstruction Using Rotation Systems},
  author={Cui, Ruiqi and G{\ae}de, Emil Toftegaard and Rotenberg, Eva and Kobbelt, Leif and B{\ae}rentzen, J Andreas},
  journal={arXiv preprint arXiv:2402.01893},
  year={2024}
}
```
