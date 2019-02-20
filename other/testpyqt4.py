from PyQt4 import QtGui  # (the example applies equally well to PySide)
import pyqtgraph as pg

## Always start by initializing Qt (only once per application)
app = QtGui.QApplication([])

## Define a top-level widget to hold everything
w = QtGui.QWidget()

## Create some widgets to be placed inside
btn = QtGui.QPushButton('press me')
text = QtGui.QLineEdit('enter text')
listw = QtGui.QListWidget()
plot = pg.PlotWidget()

## Create a grid layout to manage the widgets size and position
layout = QtGui.QGridLayout()
w.setLayout(layout)

## Add widgets to the layout in their proper positions
layout.addWidget(btn, 0, 0)   # button goes in upper-left
layout.addWidget(text, 1, 0)   # text edit goes in middle-left
layout.addWidget(listw, 2, 0)  # list widget goes in bottom-left
layout.addWidget(plot, 0, 1, 3, 1)  # plot goes on right side, spanning 3 rows

## Display the widget as a new window
w.show()

## Start the Qt event loop
app.exec_()


'''
#Exporting a file

import pyqtgraph as pg
import pyqtgraph.exporters

# generate something to export
plt = pg.plot([1,5,2,4,3])

# create an exporter instance, as an argument give it
# the item you wish to export
exporter = pg.exporters.ImageExporter(plt.plotItem)

# set export parameters if needed
exporter.parameters()['width'] = 100   # (note this also affects height parameter)

# save to file
exporter.export('fileName.png')
'''