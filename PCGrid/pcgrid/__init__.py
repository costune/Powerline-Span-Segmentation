import os
import numpy as np
import laspy
from typing import Tuple, List

from pcgrid._C import process_points

class SpanSegment:
    def __init__(self):
        self.line_las_header = None
        self.tower_las_header = None
    
    def segment(self, line_pcd_path: str, tower_pcd_path: str) -> List[np.ndarray]:
        
        if not os.path.exists(line_pcd_path):
            raise FileNotFoundError(f"{line_pcd_path} Not Found!")
        if not os.path.exists(tower_pcd_path):
            raise FileNotFoundError(f"{tower_pcd_path} Not Found!")
        
        line_las = laspy.read(line_pcd_path)
        self.line_las_header = line_las.header
        line_pcd_x = np.array(line_las.x, dtype=np.float64)
        line_pcd_y = np.array(line_las.y, dtype=np.float64)
        line_pcd_z = np.array(line_las.z, dtype=np.float64)
        line_pcd = np.stack((line_pcd_x - self.line_las_header.offsets[0],
                              line_pcd_y - self.line_las_header.offsets[1],
                              line_pcd_z), dtype=np.float32, axis=1)

        tower_las = laspy.read(tower_pcd_path)
        self.tower_las_header = tower_las.header
        tower_pcd_x = np.array(tower_las.x, dtype=np.float64)
        tower_pcd_y = np.array(tower_las.y, dtype=np.float64)
        tower_pcd_z = np.array(tower_las.z, dtype=np.float64)
        tower_pcd = np.stack((tower_pcd_x - self.tower_las_header.offsets[0],
                               tower_pcd_y - self.tower_las_header.offsets[1],
                               tower_pcd_z), dtype=np.float32, axis=1)
        
        if line_pcd.shape[0] < 100 or tower_pcd.shape[0] < 100:
            raise RuntimeError("Not enough points")
        
        try:
            span_pcd_result = process_points(line_pcd, tower_pcd)
        except Exception as e:
            raise RuntimeError("process_points failed") from e
        
        if span_pcd_result is None:
            raise RuntimeError("No result from process_points")
        
        span_pcds = []
        for pcd in span_pcd_result:
            if pcd.shape[0] <= 0:
                continue
            new_pcd = pcd.astype(np.float64)
            new_pcd[:, 0] += self.line_las_header.offsets[0]
            new_pcd[:, 1] += self.line_las_header.offsets[1]
            span_pcds.append(new_pcd)
        
        return span_pcds