import numpy as np


class MathOperationsGrid:
    def __init__(self, nx, ny, dx, dy):
        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy   

    def compute_divergence(self,vx, vy):
        div_x = np.gradient(vx, axis=0) / self.dx
        div_y = np.gradient(vy, axis=1) / self.dy
        return div_x + div_y
        
    def compute_curl(self, vx, vy):
        curl_x = np.gradient(vy, axis=0) / self.dx - np.gradient(vx, axis=1) / self.dy
        return curl_x

# Example usage:
# vx and vy are arrays representing the x and y components of the velocity field
# dx and dy are the grid spacing in the x and y directions

# Compute divergence
# divergence = compute_divergence(vx, vy, dx, dy)

# Compute curl
# curl = compute_curl(vx, vy, dx, dy)
