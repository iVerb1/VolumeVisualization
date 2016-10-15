/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package volvis;

/**
 *
 * @author s113958
 */
public class Voxel {
    
    public Voxel neighborBelow = null;
    public Voxel neighborAbove = null;
    public Voxel neighborLeft = null;
    public Voxel neighborRight = null;
    public Voxel neighborFront = null;
    public Voxel neighborBack = null;
    private short value;
    
    // coordinate of the front-top-right corner of the voxel
    private double[] origin;
    
    private static final double DIM_X = 1;
    private static final double DIM_Y = 1;
    private static final double DIM_Z = 1;
    
    public Voxel(double[] origin, short value) {
        this.origin = origin;
        this.value = value;
    }
    
    public short getValue() {
        return this.value;
    }
}
