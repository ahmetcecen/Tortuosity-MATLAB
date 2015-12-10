# Tortuosity Distribution Calculation of 3D Materials Datasets

Below should be concise enough legend of the calculation assumptions. The tortuosity is defined as the ratio of the actual (red) path to macroscale (black) path, calculated as a distribution by considering every entry point available in the entry surface. The connectivity neighborhood is assumed to be 6 nearest neighbors, as most FEM or FVM simulations will have the same connectivity.

![Legend](/Readme.png)

[TD]=Tort3D(Alpha) finds the tortuosity distribution of a properly weighted 3D dataset Alpha. Every voxel in the dataset Alpha must have a value representing the tortuosity of the phase that the voxel represents. TD is a list of calculated tortuosity values. 

[TD Tormap]=Tort3D(Alpha) will also generate a non-unique visualization of the calculated paths. 

# Pending License Approvals
