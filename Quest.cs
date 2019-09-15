using System.Collections;
using System.Collections.Generic;
using System;
using System.Numerics;
using UnityEngine;
using Random = System.Random;
/*
 * t is measured on an abstract "intensity" scale
 * tGuessSd = 3 and range = 5 for convention
 *
 * tGuess is your prior threshold estimate.
 * tGuessSd is the sd you assign to that guess. Be generous
 * pThreshold is your threshold criterion expressed as probability of response == 1.
 * beta, delta, gamma are the parameters of a Weibull psychometric function.
 * beta controls the steepness of the psychometric function. Typically 3.5
 * delta is the fraction of trials on which the observer presses blindly. Typically 0.01
 * gamma is the fraction of trials that will generate response 1 when intensity == -inf
 * grain is the quantization (step size) of the interval table. e.g. 0.01
 * range is the intensity difference between the largest and smallest intensity that the internal table can store
 **/

public class Quest {

#pragma warning disable IDE0044 // Add readonly modifier
    double tGuess;
    double tGuessSd;
    double pThreshold;
    double xThreshold;

    double beta;
    double delta;
    double gamma;

    double grain;
    double range;
    double plotIt;

    double dim;

    int[] i;
    double[] x;
    double[] pdf;

    double[] i2;
    double[] x2;
    double[] p2;

    bool updatePdf;
    bool warnPdf;
    bool normalizePdf;

    double quantileOrder;

    List<double> intensity;
    List<int> response;

    double[,] s2; // Unknown Right Now...
#pragma warning restore IDE0044 // Add readonly modifier

    public Quest(double tGuess, double tGuessSd, double pThreshold,
            double beta, double delta, double gamma,
            double grain = 0.01, double range = -1) {
#pragma warning disable RECS0117 // Local variable has the same name as a member and hides it
        double dim;
#pragma warning restore RECS0117 // Local variable has the same name as a member and hides it
        if (range <= 0)
            dim = 500;
        else {
            dim = range / grain;
            dim = 2 * Math.Ceiling(dim / 2);
        }

        this.updatePdf = true;
        this.warnPdf = true;
        this.normalizePdf = true;


        this.tGuess = tGuess;
        this.tGuessSd = tGuessSd;
        this.pThreshold = pThreshold;

        this.beta = beta;
        this.delta = delta;
        this.gamma = gamma;

        this.grain = grain;
        this.dim = dim;

        this.intensity = new List<double>();
        this.response = new List<int>();

        this.Recompute();
    }

    /*
     * Call this immediately after changing a parameter of the psychometric function.
     * QuestRecompute uses the specified parameters in "q" to recompute psychometric function.
     * It the uses the newly computed psychometric function and the history in q.intensity and
     * q.response to recompute the pdf. (QuestRecompute does nothing if q.updatePdf is false.)
     *
     * dim is the number of the distinct intensities that the internal tables in q can store.
     * The vector, of length dim, with increment size grain, will be centered on the initial guess
     * tGuess, i.e tGuess+[-range/2:grain/range/2]. QUEST assumes that intensities outside of this
     * interval have zero prior probability, i.e. they are impossible values for threshold. The
     * cost of making dim too big is some extra storage and computation, which are usually
     * negligible. Getting out-of-range warnings from QuestUpdate is one possible indication
     * that your stated range is too small.
     *
     **/

    public void Recompute() {
        if (!this.updatePdf)   return;

        if (this.gamma > this.pThreshold)
            this.gamma = 0.5;

        this.i = new int[(int)this.dim + 1];
        this.x = new double[(int)this.dim + 1];
        this.pdf = new double[(int)this.dim + 1];

        for (int cur = 0; cur < this.i.Length; cur ++)
            this.i[cur] = (int)(-this.dim / 2) + cur;

        for (int cur = 0; cur < this.x.Length; cur++)
            this.x[cur] = this.grain * this.i[cur];

        double sum = 0;
        for (int cur = 0; cur < this.pdf.Length; cur++) {
            this.pdf[cur] = Math.Exp(-0.5*Math.Pow(this.x[cur]/this.tGuessSd, 2));
            sum += this.pdf[cur];
        }

        for (int cur = 0; cur < this.pdf.Length; cur++)
            this.pdf[cur] /= sum;


        this.i2 = new double[(int)this.dim * 2 + 1];
        this.x2 = new double[(int)this.dim * 2 + 1];
        this.p2 = new double[(int)this.dim * 2 + 1];

        for (int cur = 0; cur < this.i2.Length; cur++)
            this.i2[cur] = -this.dim + cur;

        for (int cur = 0; cur < this.x2.Length; cur++)
            this.x2[cur] = this.grain * this.i2[cur];

        // Weibull Function
        // p2=delta*gamma+(1-delta)*(1-(1-gamma)*exp(-10.^(beta*(x-xThreshold))))
        for (int cur = 0; cur < this.p2.Length; cur++) {
            this.p2[cur] = this.delta * this.gamma + (1 - this.delta) * (1 - (1 - this.gamma)
                    * Math.Exp(-Math.Pow(10, this.beta * this.x2[cur])));
        }

        if (Math.Min(this.p2[0], this.p2[this.p2.Length - 1]) > this.pThreshold || Math.Max(this.p2[0], this.p2[this.p2.Length - 1]) < this.pThreshold)
            Debug.LogError("Psychometric function omits the threshold!!!");

        for (int cur = 0; cur < this.p2.Length; cur++)
            if (double.IsInfinity(this.p2[cur]))
                Debug.LogError("Psychometric function p2 is not finite");


        List<int> index = new List<int>();
        for (int cur = 0; cur < this.p2.Length - 1; cur++) {
            double diff = this.p2[cur + 1] - this.p2[cur];
            if (diff > 0)
                index.Add(cur);
        }
        if (index.Count < 2)
            Debug.LogError("Psychometric function p2 has fewer monotonic points!");

        this.xThreshold = Interpolation(GetElems(this.p2, index), GetElems(this.x2, index), this.pThreshold);

        if (double.IsInfinity(this.xThreshold))
            Debug.LogError("Psychometric function has no threshold");

        for (int cur = 0; cur < this.p2.Length; cur++)
            this.p2[cur] = this.delta * this.gamma + (1 - this.delta) * (1 - (1 - gamma)
                        * Math.Exp(-Math.Pow(10, this.beta * (this.x2[cur] + this.xThreshold))));

        for (int cur = 0; cur < this.p2.Length; cur++)
            if (double.IsInfinity(this.p2[cur]))
                Debug.LogError("Psychometric function p2 is not finite");

        this.s2 = new double[this.p2.Length, 2];
        for (int cur = 0; cur < this.p2.Length; cur++) {
            s2[cur, 0] = 1 - this.p2[this.p2.Length - 1 - cur];
            s2[cur, 1] = this.p2[this.p2.Length - 1 - cur];
        }
        // Until this part is true

        double pL = this.p2[0];
        double pH = this.p2[this.p2.Length - 1];
        double pE = pH * Math.Log(pH + double.Epsilon)
            - pL * Math.Log(pL + double.Epsilon)
            + (1 - pH + double.Epsilon) * Math.Log(1 - pH + double.Epsilon)
            - (1 - pL + double.Epsilon) * Math.Log(1 - pL + double.Epsilon);
        pE = 1 / (1 + Math.Exp(pE / (pL - pH)));
        this.quantileOrder = (pE - pL) / (pH - pL);

        for (int cur = 0; cur < this.pdf.Length; cur++)
            if (double.IsInfinity(this.pdf[cur]))
                Debug.LogError("Prior pdf is not finite!");

        double sumPdf = 0;

        for (int cur = 0; cur < this.intensity.Count; cur++) {
            double inten = Math.Max(-1e10, Math.Min(1e10, this.intensity[cur]));
            int[] ii = new int[this.i.Length];
            for (int j = 0; j < ii.Length; j++)
                ii[j] = this.pdf.Length + this.i[j]
                    - (int)(Math.Round((inten - this.tGuess) / this.grain));
            if (ii[0] < 0) {
                int sm = ii[0];
                for (int j = 0; j < ii.Length; j++)
                    ii[j] = ii[j] - sm;
            }
            if (ii[ii.Length - 1] >= this.s2.GetLength(0)) {
                int sm = ii[ii.Length - 1];
                for (int j = 0; j < ii.Length; j++)
                    ii[j] = ii[j] + this.s2.GetLength(0) - sm - 1;
            }
            sumPdf = 0;
            for (int j = 0; j < this.pdf.Length; j++) {
                this.pdf[j] = this.pdf[j] * this.s2[ii[j], this.response[cur]];
                sumPdf += this.pdf[j];
            }

            if (this.normalizePdf && cur % 100 == 0)
                for (int j = 0; j < this.pdf.Length; j++)
                    this.pdf[j] /= sumPdf;
        }

        if (this.normalizePdf) {
            sumPdf = 0;
            for (int cur = 0; cur < this.pdf.Length; cur++)
                sumPdf += this.pdf[cur];
            for (int cur = 0; cur < this.pdf.Length; cur++)
                this.pdf[cur] /= sumPdf;
        }

        for (int cur = 0; cur < this.pdf.Length; cur++)
            if (double.IsInfinity(this.pdf[cur]))
                Debug.LogError("Pdf is not finite!");
    }

    public void Update(double intensity, int response) {
        if (response < 0 || response >= this.s2.GetLength(1))
            Debug.LogError("Response out of range");

        if (this.updatePdf) {
            double sumPdf = 0;
            double inten = Math.Max(-1e10, Math.Min(1e10, intensity));
            int[] ii = new int[this.i.Length];
            for (int j = 0; j < ii.Length; j++)
                ii[j] = this.pdf.Length + this.i[j]
                    - (int)(Math.Round((inten - this.tGuess) / this.grain));

            //Reason Found
            if (ii[0] < 0) {
                int sm = ii[0];
                for (int j = 0; j < ii.Length; j++)
                    ii[j] = ii[j] - sm;
            }
            // Reason Found
            if (ii[ii.Length - 1] >= this.s2.GetLength(0)) {
                int sm = ii[ii.Length - 1];
                for (int j = 0; j < ii.Length; j++)
                    ii[j] = ii[j] + this.s2.GetLength(0) - sm - 1;
            }

            for (int j = 0; j < this.pdf.Length; j++) {
                this.pdf[j] = this.pdf[j] * this.s2[ii[j], (response + 1) % 2];
                sumPdf += this.pdf[j];
            }

            if (this.normalizePdf)
                for (int cur = 0; cur < this.pdf.Length; cur++)
                    this.pdf[cur] /= sumPdf;
        }

        this.intensity.Add(intensity);
        this.response.Add(response);
    }

    /*
     * Get the mean threshold estimate.
     **/

    public double Mean() {
        double sumVal = 0;
        double sumPdf = 0;

        for (int cur = 0; cur < this.pdf.Length; cur++) {
            sumVal += this.pdf[cur] * this.x[cur];
            sumPdf += this.pdf[cur];
        }
        return this.tGuess + sumVal / sumPdf;
    }

    public Tuple<double, double> Mode() {
        double mx = double.NegativeInfinity;
        int mxIndex = -1;
        double t;
        for (int j = 0; j < this.pdf.Length; j++) {
            if (this.pdf[j] > mx) {
                mx = this.pdf[j];
                mxIndex = j;
            }
        }
        t = this.x[mxIndex] + this.tGuess;
        return new Tuple<double, double>(t, mx);
    }

    /*
     * Gets a quantile of the pdf
     **/

    public double Quantile(double quantileOrder = -1) {
        if (quantileOrder == -1)
            quantileOrder = this.quantileOrder;

        if (quantileOrder > 1 || quantileOrder < 0)
            Debug.LogError("quantileOrder is outside range 0 to 1.");

        double[] p = new double[this.pdf.Length];
        double sum = 0;


        for (int cur = 0; cur < p.Length; cur++) {
            p[cur] = sum + this.pdf[cur];
            sum += this.pdf[cur];
        }

        if (double.IsInfinity(sum))
            Debug.LogError("Pdf is not finite");
        if (Math.Abs(sum) < Double.Epsilon)
            Debug.LogError("Pdf is all zero");

        if (quantileOrder < p[0])
            return this.tGuess + this.x[0];

        if (quantileOrder > p[p.Length - 1])
            return this.tGuess + this.x[this.x.Length - 1];

        List<int> index = new List<int>();
        double current = -1;
        for (int cur = 0; cur < p.Length; cur++) {
            double diff = p[cur] - current;
            if (diff > 0)
                index.Add(cur);
            current = p[cur];
        }
        if (index.Count < 2)
            Debug.LogError("Psychometric function p2 has fewer monotonic points!");


        double t = this.tGuess + Interpolation(GetElems(p, index), GetElems(this.x, index), p[p.Length-1] * quantileOrder);
        return t;
    }

    public int Simulate(double tTest, double tActual) {
        double t = Math.Min(this.x2[this.x2.Length - 1],
                Math.Max(tTest - tActual, this.x2[0]));
        Random random = new Random();
        double interp = Interpolation(this.x2, this.p2, t);
        int response = interp > random.NextDouble() ? 0 : 1;
        Debug.Log(t + " " + interp);
        return response;
    }


    public double Sd() {
        double sum = 0;
        double probSqSum = 0;
        double probSum = 0;

        for (int cur = 0; cur < this.pdf.Length; cur++) {
            sum += this.pdf[cur];
            probSqSum += this.pdf[cur] * Math.Pow(this.x[cur], 2);
            probSum += this.pdf[cur] * this.x[cur];
        }
        Debug.Log("probSq : " + probSqSum / sum);
        Debug.Log("prob : " + Math.Pow(probSum / sum, 2));

        return Math.Sqrt(probSqSum / sum - Math.Pow(probSum / sum, 2));
    }

#pragma warning disable RECS0082 // Parameter has the same name as a member and hides it
    private double Interpolation(double[] x, double[] v, double xq)
#pragma warning restore RECS0082 // Parameter has the same name as a member and hides it
    {
        List<KeyValuePair<double, double>> function = new List<KeyValuePair<double, double>>();
        for (int cur = 0; cur < x.Length; cur++)
            function.Add(new KeyValuePair<double, double>(x[cur], v[cur]));

        function.Sort((a,b) => (a.Value.CompareTo(b.Value)));

        for (int cur = 0; cur < function.Count; cur++)
            if (function[cur].Key > xq) {
                double xb = function[cur].Key;
                double xa = function[cur - 1].Key;
                double yb = function[cur].Value;
                double ya = function[cur - 1].Value;
                return ya + (yb - ya) * (xq - xa) / (xb - xa);
            }


        return -1;
    }

    private double[] GetElems(double[] elems, List<int> indexes) {
        double[] selectedElems = new double[indexes.Count];

        for (int cur = 0; cur < indexes.Count; cur++)
            selectedElems[cur] = elems[indexes[cur]];

        return selectedElems;
    }

    public override string ToString() {
        return "updatePdf: " + this.updatePdf + "\n" +
            "warnPdf: " + this.warnPdf + "\n" +
            "normalizePdf: " + this.normalizePdf + "\n" +
            "tGuess: " + this.tGuess + "\n" +
            "tGuessSd: " + this.tGuessSd + "\n" +
            "pThreshold: " + this.pThreshold + "\n" +
            "beta: " + this.beta + "\n" +
            "delta: " + this.delta + "\n" +
            "gamma: " + this.gamma + "\n" +
            "grain: " + this.grain + "\n" +
            "dim: " + this.dim + "\n";
    }
}

