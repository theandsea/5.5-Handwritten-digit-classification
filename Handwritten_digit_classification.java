import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Handwritten_digit_classification {
    public static void main(String[] args) throws Exception {
        x3y0 = path_int3arr("train3.txt");
        x5y1 = path_int3arr("train5.txt");
        double rate = 1.0 / 5.0 / (x3y0.length + x5y1.length);

        // training
        w = new double[64];
        double[] gradient;
        double A = 0;
        int times = 0;
        double log_lasttime = 0;
        double log_this = 0;
        ArrayList<double[]> log_record = new ArrayList<>();
        do {
            gradient = w_gradient();
            w_backward(gradient, rate);
            A = gradient_A(gradient);

            log_lasttime = log_this;
            log_this = w_forward();
            if (times % 10 == 0) {
                System.out.println("============================================");
                System.out.println(times + "___" + log_this + "____" + A);
            }
            log_record.add(new double[]{log_this});
            times++;
        } while (Math.abs(log_lasttime - log_this) > .001);

        // plot
        plot.record_draw_plot(log_record, "plot for the log-likelihood", "times", "logP");

        // the weight
        System.out.println("====================weight matrix 8*8========================");
        for (int i = 0, il = w.length; i < il; i++) {
            System.out.printf("%3.6f\t\t", w[i]);
            if ((i + 1) % 8 == 0)
                System.out.println();
        }

        // percent error rate in training data
        System.out.println("====================error rate in the training data========================");
        w_errorrate(x3y0, x5y1);

        // percent error rate in test data
        System.out.println("====================error rate in the testing data========================");
        int[][] x3y0_test = path_int3arr("test3.txt");
        int[][] x5y1_test = path_int3arr("test5.txt");
        w_errorrate(x3y0_test, x5y1_test);
    }

    public static double w_errorrate(int[][] x3_y0, int[][] x5_y1) {
        int error = 0;
        double proba_1;
        // y=0
        for (int t = 0, T = x3_y0.length; t < T; t++) {
            proba_1 = sigma(w, x3_y0[t]); // P(y=1|X)
            if (proba_1 > 0.5)
                error++;
        }
        // y=1
        for (int t = 0, T = x5_y1.length; t < T; t++) {
            proba_1 = sigma(w, x5_y1[t]);
            if (proba_1 < 0.5)
                error++;
        }

        int sumtimes = x3_y0.length + x5_y1.length;
        double errorrate = (double) error / (double) sumtimes;
        System.out.println("percent error rate___" + errorrate);
        System.out.println("error count___" + error);
        System.out.println("total count(T)___" + sumtimes);
        return errorrate;
    }

    public static double w_forward() {
        double logsum = 0;
        // y=0
        for (int t = 0, T = x3y0.length; t < T; t++) {
            logsum += Math.log(1 - sigma(w, x3y0[t]));
        }
        // y=1
        for (int t = 0, T = x5y1.length; t < T; t++) {
            logsum += Math.log(sigma(w, x5y1[t]));
        }

        return logsum;
    }

    public static void w_backward(double[] gradient, double rate) {
        for (int i = 0, il = w.length; i < il; i++) {
            w[i] += rate * gradient[i];
        }
    }

    public static double gradient_A(double[] gradient) {
        double A = 0;
        for (int i = 0, il = gradient.length; i < il; i++) {
            A += gradient[i] * gradient[i];
        }
        return A;
    }

    public static double[] w_gradient() {
        double[] gradient = new double[64];
        for (int a = 0; a < 64; a++) {
            double sum = 0;
            for (int t = 0, T = x3y0.length; t < T; t++) {
                sum += (0 - sigma(w, x3y0[t])) * x3y0[t][a];
            }
            for (int t = 0, T = x5y1.length; t < T; t++) {
                sum += (1 - sigma(w, x5y1[t])) * x5y1[t][a];
            }
            gradient[a] = sum;
        }
        return gradient;
    }

    public static double sigma(double[] weight, int[] xx) {
        // z
        double z = 0;
        for (int i = 0, il = weight.length; i < il; i++) {
            z += weight[i] * xx[i];
        }
        // sigma
        double sigma = 1.0 / (1 + Math.exp(-z));
        return sigma;
    }

    public static double[] w = null;
    public static int[][] x3y0 = null;
    public static int[][] x5y1 = null;

    public static int[][] path_int3arr(String path) throws Exception {
        String str = path_str(path);
        String[] seg = str.split("\r\n");
        int[][] x = new int[seg.length][];
        String[] int_str = null;
        for (int i = 0, il = seg.length; i < il; i++) {
            int_str = seg[i].split(" ");
            if (int_str.length == 64) {
                x[i] = new int[int_str.length];
                for (int j = 0, jl = int_str.length; j < jl; j++) {
                    x[i][j] = Integer.parseInt(int_str[j]);
                }
            } else
                throw new Exception("error size of data__" + seg[i].length());
        }
        return x;
    }

    public static String path_str(String path) throws IOException {
        StringBuffer buffer = new StringBuffer();
        BufferedReader bf = new BufferedReader(new FileReader(path));
        String s = null;
        while ((s = bf.readLine()) != null) {//使用readLine方法，一次读一行
            buffer.append(s.trim() + "\r\n");
        }

        return buffer.toString();
    }
}
