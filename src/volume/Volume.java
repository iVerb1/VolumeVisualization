/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author michel
 */
public class Volume {

    public Volume(int xd, int yd, int zd) {
        data = new short[xd * yd * zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }

    public Volume(File file) {

        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }

    }

    public short getVoxel(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }
    
    public void setVoxel(int x, int y, int z, short value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, short value) {
        data[i] = value;
    }

    public short getVoxel(int i) {
        return data[i];
    }

    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

    public short getMinimum() {
        short minimum = data[0];
        for (int i = 0; i < data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }

    public short getMaximum() {
        short maximum = data[0];
        for (int i = 0; i < data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }

    public int[] getHistogram() {
        return histogram;
    }

    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i = 0; i < data.length; i++) {
            histogram[data[i]]++;
        }
    }
    
    public boolean validVoxelPosition(int[] position) {
        return 0 <= position[0] && position[0] < getDimX() &&
                0 <= position[1] && position[1] < getDimY() &&
                0 <= position[2] && position[2] < getDimZ();
    }
    
    public boolean inBoundingBox(double[] coordinate) {
        return 0 <= coordinate[0] && coordinate[0] < getDimX() &&
                0 <= coordinate[1] && coordinate[1] < getDimY() &&
                0 <= coordinate[2] && coordinate[2] < getDimZ();
    }
    
    public int[] getVoxelNeighborPosition(int[] voxelPos, int[] direction) {
        int[] position = {voxelPos[0] + direction[0], voxelPos[1] + direction[1], voxelPos[2] + direction[2]};
        
        if (validVoxelPosition(position))
            return position;
        else
            return null;
    };
    
    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;
}
