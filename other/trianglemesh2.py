import numpy as np
from mayavi import mlab

# Create coordinate, triangular and color data
n = 29
t = np.linspace(-np.pi, np.pi, n)
z = np.exp(1j * t)
x = z.real.copy()
y = z.imag.copy()
z = np.zeros_like(x)

triangles = np.array([(0, i, i + 1) for i in range(1, n)])
x = np.r_[0, x]
y = np.r_[0, y]
z = np.r_[1, z]
face_values = [0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 0, 11, 0, 12, 0, 13, 0, 14]

# Create black wireframe mesh
wire_mesh = mlab.triangular_mesh(x, y, z, triangles, scalars=None, line_width=0.1, representation='wireframe',
                                 color=(0, 0, 0))

# Create face coloring
cell_data = wire_mesh.mlab_source.dataset.cell_data
cell_data.scalars = face_values
cell_data.scalars.name = 'Cell data'
cell_data.update()

# Plot triangular mesh with face coloring
mesh = mlab.pipeline.set_active_attribute(wire_mesh, cell_scalars='Cell data')
surf = mlab.pipeline.surface(mesh, colormap='jet')

# Retrieve the LUT colormap of the surf object. (256x4)
# this is an array of (R, G, B, A) values (each in range 0-255)
lut = surf.module_manager.scalar_lut_manager.lut.table.to_array()

# Modify alpha channel lut
lut[:, -1] = np.ones((lut.shape[0]))*255

# Add first row white color with full opacity
lut_mod = np.vstack((np.array([1, 1, 1, 0]), lut))

# Now use this colormap and show colorbar
surf.module_manager.scalar_lut_manager.lut.table = lut_mod
mlab.colorbar(surf)
mlab.show()