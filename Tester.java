import org.apache.commons.math4.legacy.optim.univariate.UnivariatePointValuePair;

import java.util.ArrayList;

public class Tester {
    int N;
    int K;
    double[] expGrowth;
    double[] actGrowth;
    ArrayList<Double> results;
    ArrayList<ArrayList<Double>> train;
    ArrayList<ArrayList<Double>> test;
    public Tester(int n, int k, ArrayList<ArrayList<Double>> trainData, ArrayList<ArrayList<Double>> testData, ArrayList<Double> outcome)// each arraylist of doubles is a single observation sequence, on each of these the algorithm is run and makes (or does not make) a trade
    {
        N = n;
        K = k;
        train = trainData;
        test = testData;
        results = outcome;
    }
    public void test()
    {
        Algorithm alg = new Algorithm(N,K,train);
        alg.parse();
        alg.train();
        Algorithm parse = new Algorithm(N,K,test,alg.m,alg.M); // just to use Algorithm's parse method
        parse.parse();
        expGrowth = new double[parse.Y.length];
        actGrowth = new double[parse.Y.length];
        for(int i = 0;i<parse.Y.length;i++)
        {
            UnivariatePointValuePair p = alg.kelly(parse.Y[i]);
            if(p.getValue()>0)
            {
                expGrowth[i] = Math.exp(p.getValue());
                actGrowth[i] = (1 - p.getPoint()) + p.getPoint() * (results.get(i) / test.get(i).get(test.get(i).size() - 1));
            }
            else {
                expGrowth[i] = 1;
                actGrowth[i] = 1;
            }
        }
    }
}
