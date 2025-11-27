echo "Starting cleaning..."
rm airfoil
rm -rf build/
mkdir build/
cd build/
echo "Cleaning completed!"
echo "Starting cmake..."
cmake ..
echo "Cmake completed!"
echo "Starting make..."
make -j
echo "make completed!"

