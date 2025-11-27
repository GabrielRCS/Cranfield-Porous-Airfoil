echo "Starting cleaning..."
rm ./airfoil;
cd ./build/
echo "Cleaning completed!"
echo "Cleaning make..."
make clean;
echo "Cleaning make completed!"
echo "Starting make..."
make -j
echo "make completed!"
