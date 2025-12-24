# Build Instructions (Summary)
- Clone CESM3_0_alpha07f.
- Copy these modified CAM files into components/atm/cam/src/physics/cam/.
- Ensure FTorch library interface is located in libraries/FTorch directory
- Add namelist entries in namelist_definition.xml.
-  Rebuild CAM and run your case.

./case.setup
./preview_namelists
./case.build --clean-all
./case.build
