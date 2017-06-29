"""The getmeta command."""

from .base import Base
from ternpy.plots import convexhull3d


class GetMeta(Base):
    """An tt command."""

    def run(self):
        tern_phases = self.options.get('FILE')
        file_list = self.options.get('INFILE')

        hull = convexhull3d.ConvexHullData(file_list, tern_phases)
        meta = convexhull3d.FindMetastable(hull)
        meta.find_all_decomposition()
