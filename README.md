# muon_tomography

Current workflow:
1. Put data in the right folders.  
2. Run make_thetas.py script from c_isect folder. It will run c_isect c++ project.  
c_isect project takes training_data or test_data file and outouts mean theta angle and confidence interval for all voxels.  
make_thetas.py runs c_isect for all files in test_data and training_data folders and creates test_thetas and train_thetas folders.  
3. Run unet.ipynb from python folder  