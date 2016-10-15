package volvis;

import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import javax.media.opengl.GL2;
import util.TFChangeListener;
import util.VectorMath;
import volume.Volume;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private RaycastRendererPanel panel;
    private TransferFunction tFunc;
    private TransferFunctionEditor tfEditor;
    private boolean print = false;
    private boolean opacityWeighting = false;
    
    public enum RENDERING_METHOD {

        MIP, COM
    }
    private RENDERING_METHOD method = RENDERING_METHOD.MIP;

    private Volume volume = null;
    private double sampleDistance = 1;
    private boolean interpolate = true;
    private double compositingEpsilon = 0.001;
    
    double[] fV = {53, 103};
    double[] aV = {0.5,0.5};
    
    private int dimX;
    private int dimY;
    private int dimZ;

    private BufferedImage image;
    private int resolutionScalingFactor = 1;
    private double[] viewMatrix = new double[4 * 4];

    private int numThreads = 8;
    private ExecutorService rayCastThreadPool;
    private List<Runnable> rayCastRunnables;
    private List<Future> rayCastTasks;

    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");

        rayCastTasks = new ArrayList();
        rayCastRunnables = new ArrayList();
        rayCastThreadPool = Executors.newCachedThreadPool();
    }

    public void setOpacityWeighting(boolean b) {
        this.opacityWeighting = b;
    }
    
    public void setMethod(RENDERING_METHOD method) {
        this.method = method;
    }

    public void setInterpolate(boolean interpolate) {
        this.interpolate = interpolate;
    }

    public void setSampleDistance(double distance) {
        this.sampleDistance = distance;
    }

    public void setResolutionScalingFactor(int factor) {
        this.resolutionScalingFactor = factor;
    }
    
    public void setNumberOfThreads(int numThreads) {
        this.numThreads = numThreads;
    }

    public void setVolume(Volume vol) {
        volume = vol;

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        panel.setTransferFunctionEditor(tfEditor);

        this.dimX = volume.getDimX();
        this.dimY = volume.getDimY();
        this.dimZ = volume.getDimZ();

        reconfigureThreadPool();
    }
    
    public void reconfigureThreadPool() {
        
        rayCastRunnables.clear();
        
        if (image != null) {
            for (int i = 0; i < numThreads; i++) {
                int intervalWidth = (int) Math.floor(((double) image.getWidth()) / numThreads);

                if (i == 0 && numThreads != 1) {
                    rayCastRunnables.add(new RayCastRunnable(0, (i + 1) * intervalWidth));
                } else if (i == numThreads - 1 && numThreads != 1) {
                    rayCastRunnables.add(new RayCastRunnable(i * intervalWidth, image.getWidth()));
                } else if (numThreads == 1) {
                    rayCastRunnables.add(new RayCastRunnable(0, image.getWidth()));
                } else {
                    rayCastRunnables.add(new RayCastRunnable(i * intervalWidth, (i + 1) * intervalWidth));
                }
            }
        }
    }

    @Override
    public void changed() {
        List<TransferFunction.ControlPoint> opwPoints = new ArrayList();
        
        for (TransferFunction.ControlPoint cp : this.tFunc.getControlPoints()) {
            if (cp.isOpacityWeightingPoint)
                opwPoints.add(cp);
        }
        
        this.aV = new double[opwPoints.size()];
        this.fV = new double[opwPoints.size()];
        
        int i = 0;
        for (TransferFunction.ControlPoint cp : opwPoints) {
            aV[i] = cp.color.a;
            fV[i] = cp.value;
            
            i++;
        }
        
        for (i = 0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }

    private int getOPWIndex0(double value) {
        if (this.aV.length < 2)
            return -1;
        
        for (int i = 0; i < this.fV.length - 1; i++) {
            if (value >= fV[i] && value <= fV[i+1])
                return i;
        }
        
        return -1;
    }
    
    public RaycastRendererPanel getPanel() {
        return panel;
    }

    // get a voxel from the volume data by nearest neighbor interpolation
    short getVoxel(double[] coord) {

        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);

        return getVoxel(x,y,z);
    }

    short getVoxel(int x, int y, int z) {
        if ((x >= 0) && (x < this.dimX) && (y >= 0) && (y < this.dimY)
                && (z >= 0) && (z < this.dimZ)) {
            return volume.getVoxel(x, y, z);
        } else {
            return 0;
        }
    }
    
    private void logVec(String name, double[] vec) {
        System.out.println(name + " (" + vec[0] + " __ " + vec[1] + " __ " + vec[2] + ")");
    }

    private void logVec(String name, int[] vec) {
        System.out.println(name + " (" + vec[0] + " __ " + vec[1] + " __ " + vec[2] + ")");
    }

    private void logVec(String name, boolean[] vec) {
        System.out.println(name + " (" + vec[0] + " __ " + vec[1] + " __ " + vec[2] + ")");
    }

    /**
     * origin must be in bounding box otherwise 0 are returned
     *
     * @param origin
     * @param direction
     * @return
     */
    public short getVoxelsAlongRay(double[] origin, double[] direction) {
        if (!volume.inBoundingBox(origin)) {
            return 0;
        } else {
            double max = 0;

            int[][] neighborDirections = new int[3][3];
            boolean[] hasNeighborDirection = new boolean[3];

            if (direction[0] == 0) {
                hasNeighborDirection[0] = false;
            } else {
                hasNeighborDirection[0] = true;

                neighborDirections[0][0] = (direction[0] > 0 ? 1 : -1);
                neighborDirections[0][1] = 0;
                neighborDirections[0][2] = 0;
            }

            if (direction[1] == 0) {
                hasNeighborDirection[1] = false;
            } else {
                hasNeighborDirection[1] = true;

                neighborDirections[1][0] = 0;
                neighborDirections[1][1] = (direction[1] > 0 ? 1 : -1);
                neighborDirections[1][2] = 0;
            }

            if (direction[2] == 0) {
                hasNeighborDirection[2] = false;
            } else {
                hasNeighborDirection[2] = true;

                neighborDirections[2][0] = 0;
                neighborDirections[2][1] = 0;
                neighborDirections[2][2] = (direction[2] > 0 ? 1 : -1);
            }

            List<Short> voxels = new ArrayList<Short>();

            int[] position = {(int) Math.round(origin[0]), (int) Math.round(origin[1]), (int) Math.round(origin[2])};
            double[] pointOfEntry = {origin[0], origin[1], origin[2]};
            double[] base = new double[3];

            do {
                voxels.add(volume.getVoxel(position[0], position[1], position[2]));

                base[0] = ((double) position[0]) - 0.5;
                base[1] = ((double) position[1]) - 0.5;
                base[2] = ((double) position[2]) - 0.5;

                double tX = -1;
                double tY = -1;
                double tZ = -1;

                if (hasNeighborDirection[0]) {
                    tX = (base[0] + (double) neighborDirections[0][0] - pointOfEntry[0]) / direction[0];
                }

                if (hasNeighborDirection[1]) {
                    tY = (base[1] + (double) neighborDirections[1][1] - pointOfEntry[1]) / direction[1];
                }

                if (hasNeighborDirection[2]) {
                    tZ = (base[2] + (double) neighborDirections[2][2] - pointOfEntry[2]) / direction[2];
                }

                double value = 0;
                double[] middle = new double[3];

                if (tX >= 0 && (tY < 0 || tX <= tY) && (tZ < 0 || tX <= tZ)) {
                    middle[0] = pointOfEntry[0] + direction[0] * 0.5 * tX;
                    middle[1] = pointOfEntry[1] + direction[1] * 0.5 * tX;
                    middle[2] = pointOfEntry[2] + direction[2] * 0.5 * tX;

                    pointOfEntry[0] = pointOfEntry[0] + direction[0] * tX;
                    pointOfEntry[1] = pointOfEntry[1] + direction[1] * tX;
                    pointOfEntry[2] = pointOfEntry[2] + direction[2] * tX;

                    position = volume.getVoxelNeighborPosition(position, neighborDirections[0]);
                } else if (tY >= 0 && (tX < 0 || tY <= tX) && (tZ < 0 || tY <= tZ)) {
                    middle[0] = pointOfEntry[0] + direction[0] * 0.5 * tY;
                    middle[1] = pointOfEntry[1] + direction[1] * 0.5 * tY;
                    middle[2] = pointOfEntry[2] + direction[2] * 0.5 * tY;

                    pointOfEntry[0] = pointOfEntry[0] + direction[0] * tY;
                    pointOfEntry[1] = pointOfEntry[1] + direction[1] * tY;
                    pointOfEntry[2] = pointOfEntry[2] + direction[2] * tY;

                    position = volume.getVoxelNeighborPosition(position, neighborDirections[1]);
                } else if (tZ >= 0 && (tX < 0 || tZ <= tX) && (tY < 0 || tZ <= tY)) {
                    middle[0] = pointOfEntry[0] + direction[0] * 0.5 * tZ;
                    middle[1] = pointOfEntry[1] + direction[1] * 0.5 * tZ;
                    middle[2] = pointOfEntry[2] + direction[2] * 0.5 * tZ;

                    pointOfEntry[0] = pointOfEntry[0] + direction[0] * tZ;
                    pointOfEntry[1] = pointOfEntry[1] + direction[1] * tZ;
                    pointOfEntry[2] = pointOfEntry[2] + direction[2] * tZ;

                    position = volume.getVoxelNeighborPosition(position, neighborDirections[2]);
                } else {
                    throw new RuntimeException("no new position found, should not be possible");
                }

                value = this.trilinearInterpolation(middle);

                if (value > max) {
                    max = value;
                }
            } while (position != null);

            return (short) max;
        }
    }

    /*
     * @pre position in bounding box of volume
     */
    private double trilinearInterpolation(double[] position) {
        //logVec("pos",position);

        if (!this.interpolate)
            return getVoxel(position);
                
        
        int floorX = (int) Math.floor(position[0]);
        int ceilX = (int) Math.ceil(position[0]);
        int floorY = (int) Math.floor(position[1]);
        int ceilY = (int) Math.ceil(position[1]);
        int floorZ = (int) Math.floor(position[2]);
        int ceilZ = (int) Math.ceil(position[2]);
        
        double v000 = getVoxel(floorX, floorY, floorZ);
        double v100 = getVoxel(ceilX, floorY, floorZ);
        double v010 = getVoxel(floorX, ceilY, floorZ);
        double v001 = getVoxel(floorX, floorY, ceilZ);
        double v110 = getVoxel(ceilX, ceilY, floorZ);
        double v101 = getVoxel(ceilX, floorY, ceilZ);
        double v011 = getVoxel(floorX, ceilY, ceilZ);
        double v111 = getVoxel(ceilX, ceilY, ceilZ);
        
        /*
        boolean bFloorX = floorX >= 0;
        boolean bCeilX = ceilX < this.dimX;
        boolean bFloorY = floorY >= 0;
        boolean bCeilY = ceilY < this.dimY;
        boolean bFloorZ = floorZ >= 0;
        boolean bCeilZ = ceilZ < this.dimZ;

        double v000 = bFloorX && bFloorY && bFloorZ ? volume.getVoxel(floorX + dimX * (floorY + dimY * floorZ)) : 0;
        double v100 = bCeilX && bFloorY && bFloorZ ? volume.getVoxel(ceilX + dimX * (floorY + dimY * floorZ)) : 0;
        double v010 = bFloorX && bCeilY && bFloorZ ? volume.getVoxel(floorX + dimX * (ceilY + dimY * floorZ)) : 0;
        double v001 = bFloorX && bFloorY && bCeilZ ? volume.getVoxel(floorX + dimX * (floorY + dimY * ceilZ)) : 0;
        double v110 = bCeilX && bCeilY && bFloorZ ? volume.getVoxel(ceilX + dimX * (ceilY + dimY * floorZ)) : 0;
        double v101 = bCeilX && bFloorY && bCeilZ ? volume.getVoxel(ceilX + dimX * (floorY + dimY * ceilZ)) : 0;
        double v011 = bFloorX && bCeilY && bCeilZ ? volume.getVoxel(floorX + dimX * (ceilY + dimY * ceilZ)) : 0;
        double v111 = bCeilX && bCeilY && bCeilZ ? volume.getVoxel(ceilX + dimX * (ceilY + dimY * ceilZ)) : 0;
        */
        
        double x = position[0] - (double) floorX;
        double y = position[1] - (double) floorY;
        double z = position[2] - (double) floorZ;

        double invX = 1 - x;
        double invY = 1 - y;
        double invZ = 1 - z;

        return v000 * invX * invY * invZ
                + v100 * x * invY * invZ
                + v010 * invX * y * invZ
                + v001 * invX * invY * z
                + v110 * x * y * invZ
                + v101 * x * invY * z
                + v011 * invX * y * z
                + v111 * x * y * z;

    }

    // will only do something useful if the origin is in the bounding box
    // assumes the direction vector is normalized (has length 1)
    private double[] maxTrilinearInterpolatedValue(double[] origin, double[] direction) {
        double[] result = {0,0};
        
        if (!volume.inBoundingBox(origin)) {
            return result;
        }

        double max = 0;
        double[] maxPos = new double[3];
        double step = 0;

        double[] position = new double[3];

        // forwards
        while (true) {
            position[0] = origin[0] + direction[0] * this.sampleDistance * step;
            position[1] = origin[1] + direction[1] * this.sampleDistance * step;
            position[2] = origin[2] + direction[2] * this.sampleDistance * step;

            if (!volume.inBoundingBox(position)) {
                break;
            }

            double interpolated = this.interpolate ? this.trilinearInterpolation(position) : getVoxel(position);

            if (interpolated > max) {
                maxPos = position;
                max = interpolated;
            }

            step++;
        }

        // backwards
        direction[0] = -direction[0];
        direction[1] = -direction[1];
        direction[2] = -direction[2];

        // step = 1 instead of step = 0 because we don't want to do the first voxel twice
        step = 1;

        while (true) {
            position[0] = origin[0] + direction[0] * this.sampleDistance * step;
            position[1] = origin[1] + direction[1] * this.sampleDistance * step;
            position[2] = origin[2] + direction[2] * this.sampleDistance * step;

            if (!volume.inBoundingBox(position)) {
                break;
            }

            double interpolated = this.interpolate ? this.trilinearInterpolation(position) : 0;

            if (interpolated > max) {
                maxPos = position;
                max = interpolated;
            }

            step++;
        }
        
        result[0] = max;
        
        if (this.opacityWeighting)
            result[1] = this.opacityWeight(maxPos);
        
        return result;
    }

    private double[] composite(double[] origin, double[] direction) {
        double[] argb = {0, 0, 0, 0};

        if (!volume.inBoundingBox(origin)) {
            return argb;
        }

        double step = 0;

        LinkedList<double[]> positions = new LinkedList();

        double[] position = new double[3];

        /*
        direction[0] = -direction[0];
        direction[1] = -direction[1];
        direction[2] = -direction[2];
        */

        // forwards
        while (true) {
            position[0] = origin[0] + direction[0] * this.sampleDistance * step;
            position[1] = origin[1] + direction[1] * this.sampleDistance * step;
            position[2] = origin[2] + direction[2] * this.sampleDistance * step;

            if (!volume.inBoundingBox(position)) {
                break;
            }

            double[] samplePosition = {position[0], position[1], position[2]};
            positions.addFirst(samplePosition);

            step++;
        }

        // backwards
        direction[0] = -direction[0];
        direction[1] = -direction[1];
        direction[2] = -direction[2];

        // step = 1 instead of step = 0 because we don't want to do the first voxel twice
        step = 1;

        while (true) {
            position[0] = origin[0] + direction[0] * this.sampleDistance * step;
            position[1] = origin[1] + direction[1] * this.sampleDistance * step;
            position[2] = origin[2] + direction[2] * this.sampleDistance * step;

            if (!volume.inBoundingBox(position)) {
                break;
            }

            double[] samplePosition = {position[0], position[1], position[2]};
            positions.addLast(samplePosition);

            step++;
        }

        argb[0] = 1;

        for (double[] pos : positions) {
            TFColor voxelColor = tFunc.getColor((int) this.trilinearInterpolation(pos));

            double a = this.opacityWeighting ? 0.01 * this.opacityWeight(pos) : voxelColor.a;//this.trilinearInterpolation(pos, voxelA);

            if (print) {
                logVec("pos", pos);
                System.out.println("a " + a);
            }
            
            argb[0] = a + (1-a) * argb[0];
            argb[1] = a * voxelColor.r + (1-a)*argb[1];//this.trilinearInterpolation(pos, voxelR);
            argb[2] = a * voxelColor.g + (1-a)*argb[2];//this.trilinearInterpolation(pos, voxelG);
            argb[3] = a * voxelColor.b + (1-a)*argb[3];//this.trilinearInterpolation(pos, voxelB);

            if (print) {
                logVec("argb", argb);
            }

        }

        return argb;
    }

    private double[] gradientVector(double[] position) {
        double[] vector = new double[3];
        
        double[] neighbor1 = {position[0], position[1], position[2]};
        double[] neighbor2 = {position[0], position[1], position[2]};
        
        neighbor1[0]++;
        neighbor2[0]--;
        vector[0] = 0.5 * (this.trilinearInterpolation(neighbor1) - this.trilinearInterpolation(neighbor2));
        neighbor1[0]--;
        neighbor2[0]++;
        
        neighbor1[1]++;
        neighbor2[1]--;
        vector[1] = 0.5 * (this.trilinearInterpolation(neighbor1) - this.trilinearInterpolation(neighbor2));
        neighbor1[1]--;
        neighbor2[1]++;
        
        neighbor1[2]++;
        neighbor2[2]--;
        vector[2] = 0.5 * (this.trilinearInterpolation(neighbor1) - this.trilinearInterpolation(neighbor2));
        
        return vector;
    }
    
    private double opacityWeight(double[] position) {
        double fx = this.trilinearInterpolation(position);
        
        int n = this.getOPWIndex0(fx);
        
        if (n == -1)
            return 0;
        
        int n1 = n + 1;
        
        return VectorMath.length(this.gradientVector(position)) * (aV[n1] * (fx - fV[n]) + aV[n] * (fV[n1] - fx)) / (fV[1] - fV[n]);
        
    }
    
    void slicer() {
        rayCastTasks.clear();

        for (int i = 0; i < rayCastRunnables.size(); i++) {
            rayCastTasks.add(rayCastThreadPool.submit(rayCastRunnables.get(i)));
        }

        try {
            for (Future f : rayCastTasks) {
                f.get();
            }
        } catch (InterruptedException ex) {
            System.out.println(ex);
        } catch (ExecutionException ex) {
            System.out.println(ex);
        }

    }

    class RayCastRunnable implements Runnable {

        final int begin;
        final int end;

        public RayCastRunnable(int begin, int end) {
            this.begin = begin;
            this.end = end;
        }

        @Override
        public void run() {
            // vector uVec and vVec define a plane through the origin, 
            // perpendicular to the view vector viewVec
            double[] viewVec = new double[3];
            double[] uVec = new double[3];
            double[] vVec = new double[3];
            VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
            VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
            VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

            // image is square
            int imageCenter = image.getWidth() / 2;

            double[] pixelCoord = new double[3];
            double[] volumeCenter = new double[3];
            VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

            if (resolutionScalingFactor == 1) {
                for (int j = 0; j < image.getHeight(); j++) {
                    for (int i = begin; i < end; i++) {

                        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                                + volumeCenter[0];// + viewVec[0] * viewVecScale;
                        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                                + volumeCenter[1];// + viewVec[1] * viewVecScale;
                        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                                + volumeCenter[2];// + viewVec[2] * viewVecScale;

                        // BufferedImage expects a pixel color packed as ARGB in an int
                        double[] argb = getARGB(pixelCoord, viewVec);
                        image.setRGB(i, j, getPixelColor(argb));
                    }
                }
            } else {
                int width = end - begin;
                int maxXIndex = end - (width % resolutionScalingFactor);
                int maxYIndex = image.getHeight() - (image.getHeight() % resolutionScalingFactor);
                double[][] argbList = new double[(int) Math.floor(((double) width) / resolutionScalingFactor)][];

                for (int j = 0; j < image.getHeight(); j++) {
                    for (int i = begin; i < end; i++) {

                        double[] argb;

                        int iMod = (i - begin) % resolutionScalingFactor;
                        int jMod = j % resolutionScalingFactor;

                        if (iMod == 0 && jMod == 0 && i < maxXIndex && j < maxYIndex) {

                            double iNew, jNew;
                            iNew = i + (((double) resolutionScalingFactor) - 1) / 2;
                            jNew = j + (((double) resolutionScalingFactor) - 1) / 2;

                            pixelCoord[0] = uVec[0] * (iNew - imageCenter) + vVec[0] * (jNew - imageCenter)
                                    + volumeCenter[0];// + viewVec[0] * viewVecScale;
                            pixelCoord[1] = uVec[1] * (iNew - imageCenter) + vVec[1] * (jNew - imageCenter)
                                    + volumeCenter[1];// + viewVec[1] * viewVecScale;
                            pixelCoord[2] = uVec[2] * (iNew - imageCenter) + vVec[2] * (jNew - imageCenter)
                                    + volumeCenter[2];// + viewVec[2] * viewVecScale;

                            argb = getARGB(pixelCoord, viewVec);

                            argbList[(i - begin) / resolutionScalingFactor] = argb;
                        } else if (i >= maxXIndex) {
                            argb = argbList[argbList.length - 1];
                        } else {
                            argb = argbList[(i - begin) / resolutionScalingFactor];
                        }

                        image.setRGB(i, j, getPixelColor(argb));
                    }
                }
            }
       }

    }
        
    private int getPixelColor(double[] argb) {
        int c_alpha = argb[0] <= 1.0 ? (int) Math.floor(argb[0] * 255) : 255;
        int c_red = argb[1] <= 1.0 ? (int) Math.floor(argb[1] * 255) : 255;
        int c_green = argb[2] <= 1.0 ? (int) Math.floor(argb[2] * 255) : 255;
        int c_blue = argb[3] <= 1.0 ? (int) Math.floor(argb[3] * 255) : 255;
        return (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
    }
    
     // returns {a,r,g,b} of the pixel
    private double[] getARGB(double[] origin, double[] direction) {
        
        double[] argb = new double[4];

        if (method == RENDERING_METHOD.MIP) {
            double[] interpol = maxTrilinearInterpolatedValue(origin, direction);
            TFColor voxelColor = tFunc.getColor((int)interpol[0]);

            argb[0] = this.opacityWeighting ? interpol[1] : voxelColor.a;
            argb[1] = voxelColor.r;
            argb[2] = voxelColor.g;
            argb[3] = voxelColor.b;
        } else if (method == RENDERING_METHOD.COM) {
            return this.composite(origin, direction);
        } 
        return argb;
    }

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL2.GL_LINE_SMOOTH);
        gl.glHint(GL2.GL_LINE_SMOOTH_HINT, GL2.GL_NICEST);
        gl.glEnable(GL2.GL_BLEND);
        gl.glBlendFunc(GL2.GL_SRC_ALPHA, GL2.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-this.dimX / 2.0, -this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(-this.dimX / 2.0, this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, -this.dimY / 2.0, this.dimZ / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-this.dimX / 2.0, -this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glVertex3d(-this.dimX / 2.0, this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, -this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(this.dimX / 2.0, -this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, -this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-this.dimX / 2.0, -this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glVertex3d(-this.dimX / 2.0, -this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(-this.dimX / 2.0, this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(-this.dimX / 2.0, this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-this.dimX / 2.0, this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glVertex3d(-this.dimX / 2.0, this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-this.dimX / 2.0, -this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glVertex3d(-this.dimX / 2.0, -this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, -this.dimY / 2.0, this.dimZ / 2.0);
        gl.glVertex3d(this.dimX / 2.0, -this.dimY / 2.0, -this.dimZ / 2.0);
        gl.glEnd();

        gl.glDisable(GL2.GL_LINE_SMOOTH);
        gl.glDisable(GL2.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {
        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        slicer();
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL2.GL_BLEND);
        gl.glBlendFunc(GL2.GL_SRC_ALPHA, GL2.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();

        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }

}
