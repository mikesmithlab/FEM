import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.widgets import RectangleSelector, EllipseSelector
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
import numpy as np
import functools
import meshUtils as mesh

def matplotlib_helper_plot(nodepts, elements):
    '''
    Provides a marked up view of nodes with node and element labels annotated
    Too slow to use on big datasets. For that you should use showmesh
    '''
    plt.plot(nodepts[:, 0], nodepts[:, 1], 'ro')

    x_element = np.append(nodepts[elements[11], 0], nodepts[elements[11], 0][0])
    y_element = np.append(nodepts[elements[11], 1], nodepts[elements[11], 1][0])
    for i in np.arange(np.shape(elements)[0]):
        x_element = np.append(nodepts[elements[i], 0], nodepts[elements[i], 0][0])
        y_element = np.append(nodepts[elements[i], 1], nodepts[elements[i], 1][0])
        plt.plot(x_element, y_element, 'b-')
    nodelabels = np.arange(np.shape(nodepts)[0])
    elementlabels = np.arange(np.shape(elements)[0])
    elementcoords = np.mean(nodepts[elements], axis=1)

    for label, x, y in zip(nodelabels, nodepts[:, 0], nodepts[:, 1]):
        plt.annotate(label, xy=(x, y), xytext=(0, 20), textcoords='offset points', ha='right', va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    for label, x, y in zip(elementlabels, elementcoords[:, 0], elementcoords[:, 1]):
        plt.annotate(label, xy=(x, y), xytext=(0, 0), textcoords='offset points', ha='right', va='bottom',
                     bbox=dict(boxstyle='round,pad=0.5', fc='red', alpha=0.5))

    plt.show()


def tri_mesh_select_plot(nodepts, elements, selected, show=True):
    fig, ax = _tri_mesh_plot(nodepts, elements)
    ax.plot(nodepts[selected, 0],nodepts[selected, 1],'ro')
    ax.set_title('Triangular grid')
    if show:
        plt.show()
    return fig, ax

def tri_mesh_value_plot(nodepts, elements, values, show=True):
    fig, ax = _tri_mesh_plot(nodepts, elements)
    nodeval = []
    for node_id, node in enumerate(nodepts):
        nodeval.append(values[mesh.find_tri_with_node(elements, node_id)[0]])

    ax.plot(nodepts[selected,0],nodepts[selected,1],'ro')
    ax.set_title('Triangular grid')
    if show:
        plt.show()
    return fig,ax

def _tri_mesh_plot(nodepts,elements):

    triang = mtri.Triangulation(nodepts[:, 0], nodepts[:, 1], elements)
    # Set up the figure
    fig, ax = plt.subplots()

    # Plot the triangulation.
    #ax.tricontourf(triang, nodepts)
    ax.triplot(triang, 'k.-', lw=0.5)
    return fig, ax




def select_pts(nodepts,elements, values=None, selection=None, shape='ellipse',add=True):
    '''
    This function allows you to interactively select a point, rectangle or ellipse.
    This info can then be fed to functions to remove nodes/elements or query their properties
    '''
    if selection is None:
        selection = np.zeros(np.shape(nodepts)[0], dtype=bool)



    def line_select_callback(eclick, erelease):
        'eclick and erelease are the press and release events'

        x1, y1= eclick.xdata, eclick.ydata

        x2,y2= erelease.xdata, erelease.ydata

        print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
        print(" The button you used were: %s %s" % (eclick.button, erelease.button))

    def toggle_selector(event):
        print(' Key pressed.')
        if event.key in ['Q', 'q'] and RS.active:
            print(' RectangleSelector deactivated.')
            RS.set_active(False)
        if event.key in ['A', 'a'] and not RS.active:
            print(' RectangleSelector activated.')
            RS.set_active(True)

    def rect_select_node(rect, selected):
        for index, node in enumerate(nodepts):
            if (node[0] > rect[0]) & (node[0] < (rect[0] + rect [2])):
                if (node[1] > rect[1]) & (node[1] < (rect[1] + rect [3])):
                    selected[index] = add
        return selected

    def ellipse_select_node(rect, selected):
        el = patches.Ellipse((rect[0],rect[1]),rect[2],rect[3])
        for index, node in enumerate(nodepts):
            if el.contains_point(node):
                selected[index] = add
        return selected

    if values is None:
        values = mesh.default_values(elements)


    fig, ax = tri_mesh_select_plot(nodepts, elements, selection, show=False)
    if shape == 'rectangle':
        RS = RectangleSelector(ax, line_select_callback,
                                               drawtype='box', useblit=True,
                                               button=[1, 3],  # don't use middle button
                                               minspanx=5, minspany=5,
                                               spancoords='pixels',
                                               interactive=True)

    elif shape == 'ellipse':
        RS = EllipseSelector(ax, line_select_callback,
                                               drawtype='box', useblit=True,
                                               button=[1, 3],  # don't use middle button
                                               minspanx=5, minspany=5,
                                               spancoords='pixels',
                                               interactive=True)

    plt.connect('key_press_event', toggle_selector)
    plt.show()

    if shape == 'rectangle':
        edge_centers = RS.edge_centers
        #[(x,y),width,height]
        rect = [edge_centers[0][0],edge_centers[1][1],edge_centers[0][2] - edge_centers[0][0],edge_centers[1][3] - edge_centers[1][1]]
    elif shape == 'ellipse':
        edge_centers = RS.edge_centers
        rect = [edge_centers[0][1], edge_centers[1][0], edge_centers[0][2] - edge_centers[0][0],
                edge_centers[1][3] - edge_centers[1][1]]

    if shape == 'rectangle':
        selection = rect_select_node(rect, selection)
    elif shape == 'ellipse':
        selection = ellipse_select_node(rect, selection)

    tri_mesh_select_plot(nodepts, elements, selection, show=True)


    return selection



if __name__ == '__main__':
    nodepts, elements,_ = mesh.load_mesh_from_file('/home/ppzmis/Documents/PythonScripts/FEM/test.h5')

    selected = select_pts(nodepts,elements, shape='ellipse',add=True)
    selected = select_pts(nodepts, elements, selection=selected, shape='rectangle',add=False)

    node_values = np.ones(np.shape(nodepts)[0])

    node_values[selected] = 2*np.ones(np.sum(selected))
