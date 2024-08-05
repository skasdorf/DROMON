ReadMe for DROMON, a C++ Boundary Element Method (BEM)/Method of Moments (MoM) Library by Jake J. Harmon (jake.harmon@ieee.org)


#####################3
Dependencies
#####################
sudo apt install build-essential cmake gfortran libgsl-dev libboost-dev

#these may not be necessary
sudo apt install libblas-dev liblapack-dev 

wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
| gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt update

sudo apt install intel-oneapi-mkl-devel

install intel oneapi and then run this in the main directory:
. /opt/intel/oneapi/setvars.sh

#this command needs to be run every time terminal is opened.  This can be automated by adding it to .bashrc file in #the home directory (need to show hidden files to see this)

#######################
Build and Install Instructions
#######################

# First, create the build directory:
mkdir DROMON_BUILD
cd DROMON_BUILD

# Then, in the build directory, we need to call CMake as follows:
# The build type may be changed to debug, release, release with debug info
cmake -DCMAKE_INSTALL_PREFIX=/DROMON/INSTALL/DIRECTORY -DCMAKE_BUILD_TYPE=Release /DROMON/LIB/DIRECTORY

###
Here this looks like:
cmake -DCMAKE_INSTALL_PREFIX=$HOME/Documents/DROMON -DCMAKE_BUILD_TYPE=Release ../


# While in the build directory, execute the following
cmake --build . --target install


# If the library is built and installed for release and debug,
# the appropriate version will be selected in a find_package call
# assuming a "config" search

# Note that the install directory should be set as an environment variable
# i.e., DROMON_DIR = /installed/directory


#######################
HOW TO BUILD/RUN THE TESTS
#######################
# First, build the library as in the instructions above
# Make sure "NUMDIFF" is installed (https://www.nongnu.org/numdiff/) for comparing output files!
# Then, create a build directory for the tests:
mkdir DROMON_TESTS
cd DROMON_TESTS

# Then, in this directory, we need to call CMake as follows (replace Release with Debug if necessary):
cmake -DCMAKE_BUILD_TYPE=Release ../DROMON/LIB/DIRECTORY/tests
make

###
Here I had to use this instead:
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$HOME/Documents/DROMON-main/cmake ../tests

###

# Finally, execute the following to actually run the tests
ctest

# The tests will first run the normal test, then, when included, tests the output files with numdiff
# Additional files and tests may be added by modifying the CMakeLists.txt file in the tests directory






