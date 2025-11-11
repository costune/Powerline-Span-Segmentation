import os
import glob
import numpy as np
import laspy
from typing import Tuple, List
from natsort import natsorted

from pcgrid._C import process_points

class SpanSegment:
    
    @staticmethod
    def segment(line_pcd_path: str, tower_pcd_path: str, output_dir: str) -> List[np.ndarray]:
        
        if not os.path.exists(line_pcd_path):
            raise FileNotFoundError(f"{line_pcd_path} Not Found!")
        if not os.path.exists(tower_pcd_path):
            raise FileNotFoundError(f"{tower_pcd_path} Not Found!")
        
        os.makedirs(output_dir, exist_ok=True)
        for filename in os.listdir(output_dir):
            file_path = os.path.join(output_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                os.rmdir(file_path)
                
        line_pcd_path = os.path.abspath(line_pcd_path)
        tower_pcd_path = os.path.abspath(tower_pcd_path)
        output_dir = os.path.abspath(output_dir)
        
        process_points(line_pcd_path, tower_pcd_path, output_dir)
        
        span_pcd_search_path = os.path.join(output_dir, "group_*.las")
        span_pcd_pathes = glob.glob(span_pcd_search_path)
        
        span_pcds = []
        for span_pcd_path in span_pcd_pathes:
            las = laspy.read(span_pcd_path)
            if las.header.point_count < 100:
                continue
            las_x = np.array(las.x)
            las_y = np.array(las.y)
            las_z = np.array(las.z)
            span_pcds.append(np.stack([las_x, las_y, las_z], axis=1))
        
        
        return span_pcds