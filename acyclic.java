
import com.sun.tools.corba.se.idl.constExpr.And;

import java.io.*;
import java.util.*;
import java.lang.Math;


/**
 * An implementation of Algorithms for testing acylicity
 * and counting acyclic orientations of undirected graphs.
 */

public class acyclic {

    public static double binom(int n, int m) {
        int i, j;
        double BC[][] = new double[126][126];
        if (n>=0)
            if (m>n||m<0) return 0.0;

            else {
                for(i=0;i<=n;i++) BC[i][0] = 1;
                for(i=1;i<=m;i++) BC[0][i] = 0;
                for(j=1;j<=m;j++) for(i=1;i<=n;i++)
                    BC[i][j] = BC[i-1][j-1] + BC[i-1][j];
            }

        return BC[n][m];
    }
    /** The method below will implement an Algorithm to
     * take a 2-dimensional square matrix representing a
     * directed graph, and determine whether that directed 
     * is acyclic. **/

    public static boolean isDAG (boolean digraph[][]) {   

        int n = digraph.length;

        boolean isThereSink = false;
        boolean[] sinks = new boolean[n];

        for (int i=0; i<n; i++) {

            boolean[] sink_temp = new boolean[n];
            int isSink = 0;

            for (int j=0; j<n; j++) {
                boolean c1 = digraph[i][j] == false;
                if (c1) {
                    sink_temp[j] = false;
                }
                else {
                    sink_temp[j] = true;
                }

            }

            for (int k=0; k<n; k++) {
                if (sink_temp[k]==true) {
                    isSink = isSink + 1;
                }
            }
            if (isSink==0) {
                sinks[i] = true;
            }
        }

        for (int l=0; l<n; l++) {
            isThereSink = isThereSink || sinks[l];
        }

        if (isThereSink) {
            return true;
        }

        else {
            return false;
        }
           /** placeholder **/
    }

    /** The method below will implement a simple algorithm
     * which takes an 2 dimensional matrix representing a 
     * undirected graph, and generates and returns a random
     * uniform orientation of that graph. **/
 
    public static boolean[][] uniformOrient(boolean[][] graph) {

        int n = graph.length;
        boolean graph2[][]=graph;

        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                double random = Math.random();
                if (i!=j) {
                    if (graph2[i][j] == true) {

                        if (random > 0.5) {
                            graph2[i][j] = true;
                            graph2[j][i] = false;
                        } else {
                            graph2[i][j] = false;
//                            graph2[j][i] = true;
                        }
                    }
                }
                else graph2[i][j] = false;
            }
        }

        return graph2;    /** placeholder **/
    }

    /** The method below will implement a simple algorithm
     * which takes a (positive) integer n, and a probability
     * p and return a 2-dimensional matrix for an undirected
     * graph G on n vertices, generated according to the
     * random model G_{n,p}.  **/   

    public static boolean[][] erdosRenyi(int n, double p) {
	    boolean[][] graph = new boolean[n][n];
	    Random generator = new Random();

        for (int i=0; i<n; i++) {
            for (int j = 0; j < n; j++) {

                if (i==j) {
                    graph[i][j] = false;
                }
                if (generator.nextDouble()<p) {
                    graph[i][j] = true;
                    graph[j][i] = false;
                }
                else {
                    graph[i][j] = false;
                    graph[j][i] = true;
                }
            }
        }

	    return graph;  /** placeholder **/
    }
	 
	
    /** The method below will implement a dynamic programming 
     * algorithm to exactly evaluate the expected number of 
     * acyclic orientations in the random graph model G_{n,p}.
     * This algorithm should be based on the Robinson-Stanley
     * recurrence, and should run in O(n^2) time if possible.
     **/

    public static double expErdosRenyi(int n, double p) {

        double expectation=0.0;

        double[] factorials = new double[n+1];
        factorials[0] = 1;
        for (int i=1; i<n+1; i++) {
            factorials[i] = i*factorials[i-1];
        }


        double expectCoeff = Math.pow((1-p),(factorials[n]/(factorials[n-2]*factorials[2])));

        double[] robStan = new double[n+1];

        double[] powers = new double[n*n+1];
        for (int pwr=0; pwr<n*n+1; pwr++) {
            powers[pwr] = Math.pow(1/(1-p),pwr);
        }

        double[] powersMinus1 = new double[n+2];
        for (int minus=0; minus<n+2; minus++) {
            powersMinus1[minus] = Math.pow(-1,minus);
        }

        robStan[0] = 1;
        robStan[1] = 1;

        for (int k=2; k<n+1; k++) {
            double temp=0.0;
            double thing=0.0;
            for (int i=1; i<k+1; i++) {
//                temp = Math.pow(-1,(i+1))*(factorials[k]/(factorials[k-i]*factorials[i]))*Math.pow(1/(1-p),(i*(k-i)))*robStan[k-i];
                temp = powersMinus1[i+1]*(factorials[k]/(factorials[k-i]*factorials[i]))*powers[i*k]/powers[i*i]*robStan[k-i];
                thing = thing + temp;
            }
            robStan[k]=thing;
        }

        expectation=expectCoeff*robStan[n];
        return expectation;  /** placeholder **/
    }

    /**
     * The following method should estimate the variance of
     * the number of AOs over k (graph, orientation) trials
     * with the graph being drawn from G_{n,p} (Erdos-Renyi) 
     * random model and the orientation being sampled 
     * uniformly at random for this graph.  
     * For calculating variance, you will probably use the
     * DP algorithm (coded up as expErdosRenyi) to get the
     * exact value of the expectation.  
     **/

    public static double estimErdosRenyi (int n, double p, int k) {

        double count = 0;
        double edges =0;




        for (int i = 0; i<k; i++) {

            boolean[][] af = erdosRenyi(n, p);

            af = uniformOrient(af);

            for (int j=0;j<af.length;j++ ) {
                for (int l = 0; l < af.length; l++) {
                    if (af[j][l]) {
                        edges = edges + 1;
                    }
                }
            }


            if (isDAG(af)) {
                count += 1;
            }
        }
        return Math.pow(2,edges/k)*count/k;  /** placeholder **/
    }

    /**
     * You should write a main method to run tests on your
     * different algorithms. 
     * 
     **/


    public static void main(String [] main) {
        System.out.println("Testing...");
        boolean[][] a = {{false,true,true,true},
                {false,false,true,false},
                {false,true,false,false},
                {false,false,false,false}};
        boolean[][] b = {{false,true,false,false},
                {false,false,true,false},
                {false,false,false,true},
                {true,false,false,false}};
        //System.out.println(isDAG(a));
        //System.out.println(uniformOrient(a));
        for (int z = 2; z<20; z++) {

            System.out.print(expErdosRenyi(z,0.7)+", ");
            System.out.println(estimErdosRenyi(z,0.7,300));
        }
    }
}
