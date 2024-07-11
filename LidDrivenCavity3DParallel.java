import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

public class LidDrivenCavity3DParallel {
    // Constants and parameters
    private static final int nx = 128;
    private static final int ny = 128;
    private static final int nz = 16;
    private static final double kinematicViscosity = 0.01;
    private static final double lx = 1.0;
    private static final double ly = 1.0;
    private static final double lz = 1.0;
    private static final double dx = lx / nx;
    private static final double dy = ly / ny;
    private static final double dz = lz / nz;

    private static final double Ut = 1.0;
    private static final double Ub = 0.0;
    private static final double Uback = 0.0;
    private static final double Ufront = 0.0;
    private static final double Vl = 0.0;
    private static final double Vr = 0.0;
    private static final double Vback = 0.0;
    private static final double Vfront = 0.0;
    private static final double Wleft = 0.0;
    private static final double Wright = 0.0;
    private static final double Wtop = 0.0;
    private static final double Wbottom = 0.0;

    private static final int MAX_ITERATIONS = 100000;
    private static final double TOLERANCE = 1e-6;
    private static final double OMEGA = 1.25;

    private static double[][][] u, ut, v, vt, w, wt, p, divut;
    private static double[][][] Ap, Ae, Aw, An, As, Atop, Ab;

    public static void main(String[] args) {

        if (args.length > 0) {
            try {
                int numProcessors = Integer.parseInt(args[0]);
                System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", String.valueOf(numProcessors));
            } catch (NumberFormatException e) {
                System.err.println("Invalid number of processors specified. Using default.");
            }
        }

        System.out.println("Running with parallelism level: " + ForkJoinPool.getCommonPoolParallelism());

        double Re = (Ut * lx) / kinematicViscosity;
        System.out.printf("Re = %.2f%n", Re);

        double dt1 = 0.5 / kinematicViscosity / (1.0 / dx / dx + 1.0 / dy / dy + 1.0 / dz / dz);
        double dt2 = 2.0 * kinematicViscosity / (Ut * Ut);
        double dt = Math.min(dt1, dt2);
        System.out.printf("dt = %.6f%n", dt);

        initializeArrays();
        initializePoissonSolver();

        // Time advance loop
        double t = 0.0;
        double tend = 1.0;
        int nsteps = 400;

        for (int step = 0; step < nsteps; step++) {
            applyBoundaryConditions();
            computeMomentum(dt);
            computeDivergence(dt);
            solvePressurePoisson(dt);
            updateVelocities(dt);

            t += dt;
            System.out.printf("Time = %.6f%n", t);
        }


        // Post-processing
        double[][][] velMag = computeVelocityMagnitude();
        double[][][] uCenter = interpolateToCenter(u, 'u');
        double[][][] vCenter = interpolateToCenter(v, 'v');
        double[][][] wCenter = interpolateToCenter(w, 'w');

        generateVTKFile("lid_driven_cavity_3d.vtk");

        System.out.println("Simulation completed. VTK file generated.");


    }

    private static void checkForNaN(double[][][] array, char arrayName) {
        for (int j = 0; j < array.length; j++) {
            for (int i = 0; i < array[j].length; i++) {
                for (int k = 0; k < array[j][i].length; k++) {
                    if (Double.isNaN(array[j][i][k])) {
                        System.out.println("NaN detected in " + arrayName + " at [" + j + "][" + i + "][" + k + "]");
                        return;
                    }
                }
            }
        }
    }


    private static class InitializeArraysTask extends RecursiveAction {
        private final int start, end;
        private static final int THRESHOLD = 10000;

        public InitializeArraysTask(int start, int end) {
            this.start = start;
            this.end = end;
        }

        @Override
        protected void compute() {
            if (end - start <= THRESHOLD) {
                computeDirectly();
            } else {
                int mid = (start + end) / 2;
                invokeAll(new InitializeArraysTask(start, mid),
                        new InitializeArraysTask(mid, end));
            }
        }

        private void computeDirectly() {
            for (int i = start; i < end; i++) {
                int j = i / ((nx + 2) * (nz + 2));
                int remainder = i % ((nx + 2) * (nz + 2));
                int k = remainder / (nx + 2);
                int l = remainder % (nx + 2);

                u[j][l][k] = 0.0;
                ut[j][l][k] = 0.0;
                v[j][l][k] = 0.0;
                vt[j][l][k] = 0.0;
                w[j][l][k] = 0.0;
                wt[j][l][k] = 0.0;
                p[j][l][k] = 0.0;
                divut[j][l][k] = 0.0;
            }
        }
    }

    private static void initializeArrays() {
        u = new double[ny + 2][nx + 2][nz + 2];
        ut = new double[ny + 2][nx + 2][nz + 2];
        v = new double[ny + 2][nx + 2][nz + 2];
        vt = new double[ny + 2][nx + 2][nz + 2];
        w = new double[ny + 2][nx + 2][nz + 2];
        wt = new double[ny + 2][nx + 2][nz + 2];
        p = new double[ny + 2][nx + 2][nz + 2];
        divut = new double[ny + 2][nx + 2][nz + 2];
        Ap = new double[ny + 2][nx + 2][nz + 2];
        Ae = new double[ny + 2][nx + 2][nz + 2];
        Aw = new double[ny + 2][nx + 2][nz + 2];
        An = new double[ny + 2][nx + 2][nz + 2];
        As = new double[ny + 2][nx + 2][nz + 2];
        Atop = new double[ny + 2][nx + 2][nz + 2];
        Ab = new double[ny + 2][nx + 2][nz + 2];

        int totalSize = (ny + 2) * (nx + 2) * (nz + 2);
        ForkJoinPool.commonPool().invoke(new InitializeArraysTask(0, totalSize));
    }


    private static class InitializePoissonSolverTask extends RecursiveAction {
        private final int start, end;
        private static final int THRESHOLD = 10000;

        public InitializePoissonSolverTask(int start, int end) {
            this.start = start;
            this.end = end;
        }

        @Override
        protected void compute() {
            if (end - start <= THRESHOLD) {
                computeDirectly();
            } else {
                int mid = (start + end) / 2;
                invokeAll(new InitializePoissonSolverTask(start, mid),
                        new InitializePoissonSolverTask(mid, end));
            }
        }

        private void computeDirectly() {
            for (int i = start; i < end; i++) {
                int j = i / ((nx + 2) * (nz + 2));
                int remainder = i % ((nx + 2) * (nz + 2));
                int k = remainder / (nx + 2);
                int l = remainder % (nx + 2);

                // Initialize coefficients for interior points
                if (j > 0 && j < ny + 1 && l > 0 && l < nx + 1 && k > 0 && k < nz + 1) {
                    Ae[j][l][k] = 1.0 / (dx * dx);
                    Aw[j][l][k] = 1.0 / (dx * dx);
                    An[j][l][k] = 1.0 / (dz * dz);
                    As[j][l][k] = 1.0 / (dz * dz);
                    Atop[j][l][k] = 1.0 / (dy * dy);
                    Ab[j][l][k] = 1.0 / (dy * dy);
                } else {
                    // Set coefficients to 0 for boundary points
                    Ae[j][l][k] = Aw[j][l][k] = An[j][l][k] = As[j][l][k] = Atop[j][l][k] = Ab[j][l][k] = 0.0;
                }

                // Adjust coefficients for boundary conditions
                if (l == 1) Aw[j][l][k] = 0.0;  // left wall
                if (l == nx) Ae[j][l][k] = 0.0; // right wall
                if (j == 1) Ab[j][l][k] = 0.0;  // bottom wall
                if (j == ny) Atop[j][l][k] = 0.0; // top wall
                if (k == 1) As[j][l][k] = 0.0;  // back wall
                if (k == nz) An[j][l][k] = 0.0; // front wall

                // Calculate Ap
                Ap[j][l][k] = -(Ae[j][l][k] + Aw[j][l][k] + An[j][l][k] + As[j][l][k] + Atop[j][l][k] + Ab[j][l][k]);
            }
        }
    }

    private static void initializePoissonSolver() {
        int totalSize = (ny + 2) * (nx + 2) * (nz + 2);
        ForkJoinPool.commonPool().invoke(new InitializePoissonSolverTask(0, totalSize));

        // Add a check to verify initialization
        for (int j = 0; j < ny + 2; j++) {
            for (int l = 0; l < nx + 2; l++) {
                for (int k = 0; k < nz + 2; k++) {
                    if (Double.isNaN(Ap[j][l][k]) || Double.isInfinite(Ap[j][l][k])) {
                        System.out.println("Invalid Ap value at [" + j + "][" + l + "][" + k + "]: " + Ap[j][l][k]);
                    }
                }
            }
        }
    }

    private static class BoundaryConditionTask extends RecursiveAction {
        private final int start, end;
        private static final int THRESHOLD = 1000;

        public BoundaryConditionTask(int start, int end) {
            this.start = start;
            this.end = end;
        }

        @Override
        protected void compute() {
            if (end - start <= THRESHOLD) {
                computeDirectly();
            } else {
                int mid = (start + end) / 2;
                invokeAll(new BoundaryConditionTask(start, mid),
                        new BoundaryConditionTask(mid, end));
            }
        }

        private void computeDirectly() {
            for (int i = start; i < end; i++) {
                // Apply u boundary conditions
                for (int k = 0; k < nz + 2; k++) {
                    u[i][1][k] = 0.0;
                    u[i][nx + 1][k] = 0.0;
                    u[i][0][k] = 2.0 * Ub - u[i][1][k];
                    u[i][nx + 1][k] = 2.0 * Ut - u[i][nx][k];
                }

                // Apply v boundary conditions
                for (int k = 0; k < nz + 2; k++) {
                    v[i][0][k] = 2.0 * Vl - v[i][1][k];
                    v[i][nx + 1][k] = 2.0 * Vr - v[i][nx][k];
                }
                if (i == 1 || i == ny + 1) {
                    for (int j = 0; j < nx + 2; j++) {
                        v[i][j][0] = 0.0;
                    }
                }

                // Apply w boundary conditions
                for (int j = 0; j < nx + 2; j++) {
                    w[i][j][0] = 2.0 * Wleft - w[i][j][1];
                    w[i][j][nz + 1] = 2.0 * Wright - w[i][j][nz];
                }
                if (i == 0 || i == ny + 1) {
                    for (int j = 0; j < nx + 2; j++) {
                        w[i][j][0] = 2.0 * Wbottom - w[1][j][0];
                        w[i][j][nz + 1] = 2.0 * Wtop - w[ny][j][nz + 1];
                    }
                }
            }
        }
    }


    private static void applyBoundaryConditions() {
        // Apply boundary conditions for u
        for (int j = 0; j < ny + 2; j++) {
            for (int k = 0; k < nz + 2; k++) {
                // left wall
                u[j][1][k] = 0.0;
                // right wall
                u[j][nx + 1][k] = 0.0;
            }
        }

        for (int i = 0; i < nx + 2; i++) {
            for (int k = 0; k < nz + 2; k++) {
                // top wall
                u[ny + 1][i][k] = 2.0 * Ut - u[ny][i][k];
                // bottom wall
                u[0][i][k] = 2.0 * Ub - u[1][i][k];
            }
        }

        for (int j = 0; j < ny + 2; j++) {
            for (int i = 0; i < nx + 2; i++) {
                // back wall
                u[j][i][0] = 2.0 * Uback - u[j][i][1];
                // front wall
                u[j][i][nz + 1] = 2.0 * Ufront - u[j][i][nz];
            }
        }

        // Apply boundary conditions for v
        for (int j = 0; j < ny + 2; j++) {
            for (int k = 0; k < nz + 2; k++) {
                // left wall
                v[j][0][k] = 2.0 * Vl - v[j][1][k];
                // right wall
                v[j][nx + 1][k] = 2.0 * Vr - v[j][nx][k];
            }
        }

        for (int i = 0; i < nx + 2; i++) {
            for (int k = 0; k < nz + 2; k++) {
                // top wall
                v[ny + 1][i][k] = 0.0;
                // bottom wall
                v[1][i][k] = 0.0;
            }
        }

        for (int j = 0; j < ny + 2; j++) {
            for (int i = 0; i < nx + 2; i++) {
                // back wall
                v[j][i][0] = 2.0 * Vback - v[j][i][1];
                // front wall
                v[j][i][nz + 1] = 2.0 * Vfront - v[j][i][nz];
            }
        }

        // Apply boundary conditions for w
        for (int j = 0; j < ny + 2; j++) {
            for (int k = 0; k < nz + 2; k++) {
                // left wall
                w[j][0][k] = 2.0 * Wleft - w[j][1][k];
                // right wall
                w[j][nx + 1][k] = 2.0 * Wright - w[j][nx][k];
            }
        }

        for (int i = 0; i < nx + 2; i++) {
            for (int k = 0; k < nz + 2; k++) {
                // top wall
                w[ny + 1][i][k] = 2.0 * Wtop - w[ny][i][k];
                // bottom wall
                w[0][i][k] = 2.0 * Wbottom - w[1][i][k];
            }
        }

        for (int j = 0; j < ny + 2; j++) {
            for (int i = 0; i < nx + 2; i++) {
                // back wall
                w[j][i][1] = 0.0;
                // front wall
                w[j][i][nz + 1] = 0.0;
            }
        }
    }

    private static class XMomentumTask extends RecursiveAction {
        private final int startJ, endJ, startI, endI, startK, endK;
        private final double dt;
        private static final int THRESHOLD = 1000; // Adjust this value based on your system

        public XMomentumTask(int startJ, int endJ, int startI, int endI, int startK, int endK, double dt) {
            this.startJ = startJ;
            this.endJ = endJ;
            this.startI = startI;
            this.endI = endI;
            this.startK = startK;
            this.endK = endK;
            this.dt = dt;
        }

        @Override
        protected void compute() {
            if ((endJ - startJ) * (endI - startI) * (endK - startK) <= THRESHOLD) {
                computeDirectly();
            } else {
                int midJ = (startJ + endJ) / 2;
                int midI = (startI + endI) / 2;
                int midK = (startK + endK) / 2;

                invokeAll(
                        new XMomentumTask(startJ, midJ, startI, midI, startK, midK, dt),
                        new XMomentumTask(startJ, midJ, startI, midI, midK, endK, dt),
                        new XMomentumTask(startJ, midJ, midI, endI, startK, midK, dt),
                        new XMomentumTask(startJ, midJ, midI, endI, midK, endK, dt),
                        new XMomentumTask(midJ, endJ, startI, midI, startK, midK, dt),
                        new XMomentumTask(midJ, endJ, startI, midI, midK, endK, dt),
                        new XMomentumTask(midJ, endJ, midI, endI, startK, midK, dt),
                        new XMomentumTask(midJ, endJ, midI, endI, midK, endK, dt)
                );
            }
        }

        private void computeDirectly() {
            for (int j = startJ; j < endJ; j++) {
                for (int i = Math.max(startI, 2); i < endI; i++) {
                    for (int k = startK; k < endK; k++) {
                        double ue = 0.5 * (u[j][i + 1][k] + u[j][i][k]);
                        double uw = 0.5 * (u[j][i - 1][k] + u[j][i][k]);
                        double utop = 0.5 * (u[j + 1][i][k] + u[j][i][k]);
                        double ub = 0.5 * (u[j - 1][i][k] + u[j][i][k]);
                        double vtop = 0.5 * (v[j + 1][i][k] + v[j + 1][i - 1][k]);
                        double vb = 0.5 * (v[j][i][k] + v[j][i - 1][k]);
                        double un = 0.5 * (u[j][i][k + 1] + u[j][i][k]);
                        double wn = 0.5 * (w[j][i][k + 1] + w[j][i - 1][k + 1]);
                        double us = 0.5 * (u[j][i][k - 1] + u[j][i][k]);
                        double ws = 0.5 * (w[j][i][k] + w[j][i - 1][k]);

                        double convection = -(ue * ue - uw * uw) / dx - (utop * vtop - ub * vb) / dy - (un * wn - us * ws) / dz;
                        double diffusion = kinematicViscosity * (
                                (u[j][i + 1][k] - 2.0 * u[j][i][k] + u[j][i - 1][k]) / (dx * dx) +
                                        (u[j + 1][i][k] - 2.0 * u[j][i][k] + u[j - 1][i][k]) / (dy * dy) +
                                        (u[j][i][k + 1] - 2.0 * u[j][i][k] + u[j][i][k - 1]) / (dz * dz)
                        );

                        ut[j][i][k] = u[j][i][k] + dt * (convection + diffusion);
                    }
                }
            }
        }
    }

    private static class YMomentumTask extends RecursiveAction {
        private final int startJ, endJ, startI, endI, startK, endK;
        private final double dt;
        private static final int THRESHOLD = 1000; // Adjust this value based on your system

        public YMomentumTask(int startJ, int endJ, int startI, int endI, int startK, int endK, double dt) {
            this.startJ = startJ;
            this.endJ = endJ;
            this.startI = startI;
            this.endI = endI;
            this.startK = startK;
            this.endK = endK;
            this.dt = dt;
        }

        @Override
        protected void compute() {
            if ((endJ - startJ) * (endI - startI) * (endK - startK) <= THRESHOLD) {
                computeDirectly();
            } else {
                int midJ = (startJ + endJ) / 2;
                int midI = (startI + endI) / 2;
                int midK = (startK + endK) / 2;

                invokeAll(
                        new YMomentumTask(startJ, midJ, startI, midI, startK, midK, dt),
                        new YMomentumTask(startJ, midJ, startI, midI, midK, endK, dt),
                        new YMomentumTask(startJ, midJ, midI, endI, startK, midK, dt),
                        new YMomentumTask(startJ, midJ, midI, endI, midK, endK, dt),
                        new YMomentumTask(midJ, endJ, startI, midI, startK, midK, dt),
                        new YMomentumTask(midJ, endJ, startI, midI, midK, endK, dt),
                        new YMomentumTask(midJ, endJ, midI, endI, startK, midK, dt),
                        new YMomentumTask(midJ, endJ, midI, endI, midK, endK, dt)
                );
            }
        }

        private void computeDirectly() {
            for (int j = Math.max(startJ, 2); j < endJ; j++) {
                for (int i = startI; i < endI; i++) {
                    for (int k = startK; k < endK; k++) {
                        double ve = 0.5 * (v[j][i + 1][k] + v[j][i][k]);
                        double vw = 0.5 * (v[j][i - 1][k] + v[j][i][k]);
                        double vn = 0.5 * (v[j][i][k + 1] + v[j][i][k]);
                        double vs = 0.5 * (v[j][i][k - 1] + v[j][i][k]);
                        double vtop = 0.5 * (v[j + 1][i][k] + v[j][i][k]);
                        double vb = 0.5 * (v[j - 1][i][k] + v[j][i][k]);
                        double ue = 0.5 * (u[j][i + 1][k] + u[j - 1][i + 1][k]);
                        double uw = 0.5 * (u[j][i][k] + u[j - 1][i][k]);
                        double wn = 0.5 * (w[j][i][k + 1] + w[j - 1][i][k + 1]);
                        double ws = 0.5 * (w[j][i][k] + w[j - 1][i][k]);

                        double convection = -(ue * ve - uw * vw) / dx - (vtop * vtop - vb * vb) / dy - (vn * wn - vs * ws) / dz;
                        double diffusion = kinematicViscosity * (
                                (v[j][i + 1][k] - 2.0 * v[j][i][k] + v[j][i - 1][k]) / (dx * dx) +
                                        (v[j + 1][i][k] - 2.0 * v[j][i][k] + v[j - 1][i][k]) / (dy * dy) +
                                        (v[j][i][k + 1] - 2.0 * v[j][i][k] + v[j][i][k - 1]) / (dz * dz)
                        );

                        vt[j][i][k] = v[j][i][k] + dt * (convection + diffusion);
                    }
                }
            }
        }
    }

    private static class ZMomentumTask extends RecursiveAction {
        private final int startJ, endJ, startI, endI, startK, endK;
        private final double dt;
        private static final int THRESHOLD = 1000; // Adjust this value based on your system

        public ZMomentumTask(int startJ, int endJ, int startI, int endI, int startK, int endK, double dt) {
            this.startJ = startJ;
            this.endJ = endJ;
            this.startI = startI;
            this.endI = endI;
            this.startK = startK;
            this.endK = endK;
            this.dt = dt;
        }

        @Override
        protected void compute() {
            if ((endJ - startJ) * (endI - startI) * (endK - startK) <= THRESHOLD) {
                computeDirectly();
            } else {
                int midJ = (startJ + endJ) / 2;
                int midI = (startI + endI) / 2;
                int midK = (startK + endK) / 2;

                invokeAll(
                        new ZMomentumTask(startJ, midJ, startI, midI, startK, midK, dt),
                        new ZMomentumTask(startJ, midJ, startI, midI, midK, endK, dt),
                        new ZMomentumTask(startJ, midJ, midI, endI, startK, midK, dt),
                        new ZMomentumTask(startJ, midJ, midI, endI, midK, endK, dt),
                        new ZMomentumTask(midJ, endJ, startI, midI, startK, midK, dt),
                        new ZMomentumTask(midJ, endJ, startI, midI, midK, endK, dt),
                        new ZMomentumTask(midJ, endJ, midI, endI, startK, midK, dt),
                        new ZMomentumTask(midJ, endJ, midI, endI, midK, endK, dt)
                );
            }
        }

        private void computeDirectly() {
            for (int j = startJ; j < endJ; j++) {
                for (int i = startI; i < endI; i++) {
                    for (int k = Math.max(startK, 2); k < endK; k++) {
                        double we = 0.5 * (w[j][i + 1][k] + w[j][i][k]);
                        double ww = 0.5 * (w[j][i - 1][k] + w[j][i][k]);
                        double wn = 0.5 * (w[j][i][k + 1] + w[j][i][k]);
                        double ws = 0.5 * (w[j][i][k - 1] + w[j][i][k]);
                        double wtop = 0.5 * (w[j + 1][i][k] + w[j][i][k]);
                        double wb = 0.5 * (w[j - 1][i][k] + w[j][i][k]);
                        double ue = 0.5 * (u[j][i + 1][k] + u[j][i + 1][k - 1]);
                        double uw = 0.5 * (u[j][i][k] + u[j][i][k - 1]);
                        double vn = 0.5 * (v[j + 1][i][k] + v[j + 1][i][k - 1]);
                        double vs = 0.5 * (v[j][i][k] + v[j][i][k - 1]);

                        double convection = -(we * ue - ww * uw) / dx - (wtop * vn - wb * vs) / dy - (wn * wn - ws * ws) / dz;
                        double diffusion = kinematicViscosity * (
                                (w[j][i + 1][k] - 2.0 * w[j][i][k] + w[j][i - 1][k]) / (dx * dx) +
                                        (w[j + 1][i][k] - 2.0 * w[j][i][k] + w[j - 1][i][k]) / (dy * dy) +
                                        (w[j][i][k + 1] - 2.0 * w[j][i][k] + w[j][i][k - 1]) / (dz * dz)
                        );

                        wt[j][i][k] = w[j][i][k] + dt * (convection + diffusion);
                    }
                }
            }
        }
    }


    private static void computeMomentum(double dt) {
        computeXMomentum(dt);
        computeYMomentum(dt);
        computeZMomentum(dt);
    }

    private static void computeXMomentum(double dt) {
        ForkJoinPool.commonPool().invoke(new XMomentumTask(1, ny + 1, 2, nx + 1, 1, nz + 1, dt));
    }

    private static void computeYMomentum(double dt) {
        ForkJoinPool.commonPool().invoke(new YMomentumTask(2, ny + 1, 1, nx + 1, 1, nz + 1, dt));
    }

    private static void computeZMomentum(double dt) {
        ForkJoinPool.commonPool().invoke(new ZMomentumTask(1, ny + 1, 1, nx + 1, 2, nz + 1, dt));
    }

    private static class DivergenceTask extends RecursiveAction {
        private final int startJ, endJ, startI, endI, startK, endK;
        private static final int THRESHOLD = 1000;

        public DivergenceTask(int startJ, int endJ, int startI, int endI, int startK, int endK) {
            this.startJ = startJ;
            this.endJ = endJ;
            this.startI = startI;
            this.endI = endI;
            this.startK = startK;
            this.endK = endK;
        }

        @Override
        protected void compute() {
            if ((endJ - startJ) * (endI - startI) * (endK - startK) <= THRESHOLD) {
                computeDirectly();
            } else {
                int midJ = (startJ + endJ) / 2;
                int midI = (startI + endI) / 2;
                int midK = (startK + endK) / 2;

                invokeAll(
                        new DivergenceTask(startJ, midJ, startI, midI, startK, midK),
                        new DivergenceTask(startJ, midJ, startI, midI, midK, endK),
                        new DivergenceTask(startJ, midJ, midI, endI, startK, midK),
                        new DivergenceTask(startJ, midJ, midI, endI, midK, endK),
                        new DivergenceTask(midJ, endJ, startI, midI, startK, midK),
                        new DivergenceTask(midJ, endJ, startI, midI, midK, endK),
                        new DivergenceTask(midJ, endJ, midI, endI, startK, midK),
                        new DivergenceTask(midJ, endJ, midI, endI, midK, endK)
                );
            }
        }

        private void computeDirectly() {
            for (int j = startJ; j < endJ; j++) {
                for (int i = startI; i < endI; i++) {
                    for (int k = startK; k < endK; k++) {
                        divut[j][i][k] = (ut[j][i + 1][k] - ut[j][i][k]) / dx +
                                (vt[j + 1][i][k] - vt[j][i][k]) / dy +
                                (wt[j][i][k + 1] - wt[j][i][k]) / dz;
                    }
                }
            }
        }
    }

    private static void computeDivergence(double dt) {
        ForkJoinPool.commonPool().invoke(new DivergenceTask(1, ny + 1, 1, nx + 1, 1, nz + 1));
    }

    private static void solvePressurePoisson(double dt) {
        double error = 1.0;
        int iteration = 0;

        // Calculate the right-hand side of the pressure Poisson equation
        double[][][] S = new double[ny + 2][nx + 2][nz + 2];
        for (int j = 1; j < ny + 1; j++) {
            for (int i = 1; i < nx + 1; i++) {
                for (int k = 1; k < nz + 1; k++) {
                    S[j][i][k] = divut[j][i][k] / dt;
                }
            }
        }

        while (error > TOLERANCE && iteration < MAX_ITERATIONS) {
            error = 0.0;

            for (int j = 1; j < ny + 1; j++) {
                for (int i = 1; i < nx + 1; i++) {
                    for (int k = 1; k < nz + 1; k++) {
                        double pOld = p[j][i][k];

                        p[j][i][k] = (1 - OMEGA) * p[j][i][k] + OMEGA * (
                                S[j][i][k] - (
                                        Ae[j][i][k] * p[j][i + 1][k] +
                                                Aw[j][i][k] * p[j][i - 1][k] +
                                                Atop[j][i][k] * p[j + 1][i][k] +
                                                Ab[j][i][k] * p[j - 1][i][k] +
                                                An[j][i][k] * p[j][i][k + 1] +
                                                As[j][i][k] * p[j][i][k - 1]
                                )
                        ) / Ap[j][i][k];

                        error += Math.pow(p[j][i][k] - pOld, 2);
                    }
                }
            }

            error = Math.sqrt(error / (nx * ny * nz));
            iteration++;
        }

        System.out.printf("Pressure Poisson solved in %d iterations, final error: %.6e%n", iteration, error);
    }

    private static class UpdateVelocitiesTask extends RecursiveAction {
        private final int startJ, endJ, startI, endI, startK, endK;
        private final double dt;
        private static final int THRESHOLD = 1000;

        public UpdateVelocitiesTask(int startJ, int endJ, int startI, int endI, int startK, int endK, double dt) {
            this.startJ = startJ;
            this.endJ = endJ;
            this.startI = startI;
            this.endI = endI;
            this.startK = startK;
            this.endK = endK;
            this.dt = dt;
        }

        @Override
        protected void compute() {
            if ((endJ - startJ) * (endI - startI) * (endK - startK) <= THRESHOLD) {
                computeDirectly();
            } else {
                int midJ = (startJ + endJ) / 2;
                int midI = (startI + endI) / 2;
                int midK = (startK + endK) / 2;

                invokeAll(
                        new UpdateVelocitiesTask(startJ, midJ, startI, midI, startK, midK, dt),
                        new UpdateVelocitiesTask(startJ, midJ, startI, midI, midK, endK, dt),
                        new UpdateVelocitiesTask(startJ, midJ, midI, endI, startK, midK, dt),
                        new UpdateVelocitiesTask(startJ, midJ, midI, endI, midK, endK, dt),
                        new UpdateVelocitiesTask(midJ, endJ, startI, midI, startK, midK, dt),
                        new UpdateVelocitiesTask(midJ, endJ, startI, midI, midK, endK, dt),
                        new UpdateVelocitiesTask(midJ, endJ, midI, endI, startK, midK, dt),
                        new UpdateVelocitiesTask(midJ, endJ, midI, endI, midK, endK, dt)
                );
            }
        }

        private void computeDirectly() {
            // Update u velocity
            for (int j = startJ; j < endJ; j++) {
                for (int i = Math.max(startI, 2); i < endI; i++) {
                    for (int k = startK; k < endK; k++) {
                        u[j][i][k] = ut[j][i][k] - dt * (p[j][i][k] - p[j][i - 1][k]) / dx;
                    }
                }
            }

            // Update v velocity
            for (int j = Math.max(startJ, 2); j < endJ; j++) {
                for (int i = startI; i < endI; i++) {
                    for (int k = startK; k < endK; k++) {
                        v[j][i][k] = vt[j][i][k] - dt * (p[j][i][k] - p[j - 1][i][k]) / dy;
                    }
                }
            }

            // Update w velocity
            for (int j = startJ; j < endJ; j++) {
                for (int i = startI; i < endI; i++) {
                    for (int k = Math.max(startK, 2); k < endK; k++) {
                        w[j][i][k] = wt[j][i][k] - dt * (p[j][i][k] - p[j][i][k - 1]) / dz;
                    }
                }
            }
        }
    }

    private static void updateVelocities(double dt) {
        ForkJoinPool.commonPool().invoke(new UpdateVelocitiesTask(1, ny + 1, 1, nx + 1, 1, nz + 1, dt));
    }


    private static class VelocityMagnitudeTask extends RecursiveAction {
        private final int start, end;
        private final double[][][] velMag;
        private static final int THRESHOLD = 10000;

        public VelocityMagnitudeTask(int start, int end, double[][][] velMag) {
            this.start = start;
            this.end = end;
            this.velMag = velMag;
        }

        @Override
        protected void compute() {
            if (end - start <= THRESHOLD) {
                computeDirectly();
            } else {
                int mid = (start + end) / 2;
                invokeAll(new VelocityMagnitudeTask(start, mid, velMag),
                        new VelocityMagnitudeTask(mid, end, velMag));
            }
        }

        private void computeDirectly() {
            for (int i = start; i < end; i++) {
                int j = i / ((nx + 2) * (nz + 2));
                int remainder = i % ((nx + 2) * (nz + 2));
                int k = remainder / (nx + 2);
                int l = remainder % (nx + 2);

                velMag[j][l][k] = Math.sqrt(u[j][l][k] * u[j][l][k] +
                        v[j][l][k] * v[j][l][k] +
                        w[j][l][k] * w[j][l][k]);
            }
        }
    }

    private static double[][][] computeVelocityMagnitude() {
        double[][][] velMag = new double[ny + 2][nx + 2][nz + 2];
        int totalSize = (ny + 2) * (nx + 2) * (nz + 2);
        ForkJoinPool.commonPool().invoke(new VelocityMagnitudeTask(0, totalSize, velMag));
        return velMag;
    }

    private static class InterpolateToCenterTask extends RecursiveAction {
        private final int start, end;
        private final double[][][] field, centered;
        private final char component;
        private static final int THRESHOLD = 10000;

        public InterpolateToCenterTask(int start, int end, double[][][] field, double[][][] centered, char component) {
            this.start = start;
            this.end = end;
            this.field = field;
            this.centered = centered;
            this.component = component;
        }

        @Override
        protected void compute() {
            if (end - start <= THRESHOLD) {
                computeDirectly();
            } else {
                int mid = (start + end) / 2;
                invokeAll(new InterpolateToCenterTask(start, mid, field, centered, component),
                        new InterpolateToCenterTask(mid, end, field, centered, component));
            }
        }

        private void computeDirectly() {
            for (int i = start; i < end; i++) {
                int j = i / ((nx + 1) * (nz + 1));
                int remainder = i % ((nx + 1) * (nz + 1));
                int k = remainder / (nx + 1);
                int l = remainder % (nx + 1);

                switch (component) {
                    case 'u':
                        centered[j][l][k] = 0.5 * (field[j][l][k] + field[j][l + 1][k]);
                        break;
                    case 'v':
                        centered[j][l][k] = 0.5 * (field[j][l][k] + field[j + 1][l][k]);
                        break;
                    case 'w':
                        centered[j][l][k] = 0.5 * (field[j][l][k] + field[j][l][k + 1]);
                        break;
                    default:
                        centered[j][l][k] = field[j][l][k]; // For pressure or other scalar fields
                }
            }
        }
    }

    private static double[][][] interpolateToCenter(double[][][] field, char component) {
        double[][][] centered = new double[ny + 1][nx + 1][nz + 1];
        int totalSize = (ny + 1) * (nx + 1) * (nz + 1);
        ForkJoinPool.commonPool().invoke(new InterpolateToCenterTask(0, totalSize, field, centered, component));
        return centered;
    }

    private static class WriteVTKDataTask extends RecursiveAction {
        private final int start, end;
        private final double[][][] uCenter, vCenter, wCenter, pCenter;
        private final PrintWriter writer;
        private static final int THRESHOLD = 10000;

        public WriteVTKDataTask(int start, int end, double[][][] uCenter, double[][][] vCenter,
                                double[][][] wCenter, double[][][] pCenter, PrintWriter writer) {
            this.start = start;
            this.end = end;
            this.uCenter = uCenter;
            this.vCenter = vCenter;
            this.wCenter = wCenter;
            this.pCenter = pCenter;
            this.writer = writer;
        }

        @Override
        protected void compute() {
            if (end - start <= THRESHOLD) {
                computeDirectly();
            } else {
                int mid = (start + end) / 2;
                invokeAll(new WriteVTKDataTask(start, mid, uCenter, vCenter, wCenter, pCenter, writer),
                        new WriteVTKDataTask(mid, end, uCenter, vCenter, wCenter, pCenter, writer));
            }
        }

        private void computeDirectly() {
            StringBuilder velocityData = new StringBuilder();
            StringBuilder pressureData = new StringBuilder();

            for (int i = start; i < end; i++) {
                int k = i / ((ny + 1) * (nx + 1));
                int remainder = i % ((ny + 1) * (nx + 1));
                int j = remainder / (nx + 1);
                int l = remainder % (nx + 1);

                velocityData.append(String.format("%.6e %.6e %.6e\n",
                        uCenter[j][l][k], vCenter[j][l][k], wCenter[j][l][k]));
                pressureData.append(String.format("%.6e\n", pCenter[j][l][k]));
            }

            synchronized (writer) {
                writer.print(velocityData);
                writer.print(pressureData);
            }
        }
    }


    private static void generateVTKFile(String filename) {
        try (PrintWriter writer = new PrintWriter(new FileWriter(filename))) {
            // VTK file header
            writer.println("# vtk DataFile Version 3.0");
            writer.println("Lid Driven Cavity Flow");
            writer.println("ASCII");
            writer.println("DATASET RECTILINEAR_GRID");
            writer.println("DIMENSIONS " + (nx + 1) + " " + (ny + 1) + " " + (nz + 1));

            // X coordinates
            writer.println("X_COORDINATES " + (nx + 1) + " float");
            for (int i = 0; i < nx + 1; i++) {
                writer.println(i * dx);
            }

            // Y coordinates
            writer.println("Y_COORDINATES " + (ny + 1) + " float");
            for (int j = 0; j < ny + 1; j++) {
                writer.println(j * dy);
            }

            // Z coordinates
            writer.println("Z_COORDINATES " + (nz + 1) + " float");
            for (int k = 0; k < nz + 1; k++) {
                writer.println(k * dz);
            }

            // Point data
            writer.println("POINT_DATA " + ((nx + 1) * (ny + 1) * (nz + 1)));

            // Velocity vector
            writer.println("VECTORS velocity float");
            double[][][] uCenter = interpolateToCenter(u, 'u');
            double[][][] vCenter = interpolateToCenter(v, 'v');
            double[][][] wCenter = interpolateToCenter(w, 'w');
            for (int k = 0; k < nz + 1; k++) {
                for (int j = 0; j < ny + 1; j++) {
                    for (int i = 0; i < nx + 1; i++) {
                        writer.println(uCenter[j][i][k] + " " + vCenter[j][i][k] + " " + wCenter[j][i][k]);
                    }
                }
            }

            // Pressure scalar
            writer.println("SCALARS pressure float");
            writer.println("LOOKUP_TABLE default");
            double[][][] pCenter = interpolateToCenter(p, 'p');
            for (int k = 0; k < nz + 1; k++) {
                for (int j = 0; j < ny + 1; j++) {
                    for (int i = 0; i < nx + 1; i++) {
                        writer.println(pCenter[j][i][k]);
                    }
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}