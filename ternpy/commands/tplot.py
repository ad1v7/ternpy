"""The tplot command."""

from .base import Base
from ternpy.plots import convexhull3d


class TPlot(Base):
    """An tt command."""

    def run(self):
        tern_phases = self.options.get('FILE')
        file_list = self.options.get('INFILE')

        hull = convexhull3d.ConvexHullData(file_list, tern_phases)
        convexhull3d.PlotConvexHull(hull, steps=12).plot_all()
