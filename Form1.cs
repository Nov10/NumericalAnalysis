using System.Diagnostics;

namespace WinFormsApp1
{
    public partial class Form1 : Form
    {
        //public Form1()
        //{
        //    InitializeComponent();
        //}
        public Form1()
        {
            InitializeComponent();
            h = (rangeEnd - rangeStart) / N;

            D();
            //chart1.Series.Add("");
            //for (int i = 0; i < N + 1; i++)
            //{
            //    chart1.Series[1].Points.AddXY(rangeStart + (i - 1) * h + 0.00000000001, function(rangeStart + i * h));
            //    chart1.Series[1].Points.AddXY(rangeStart + i * h, function(rangeStart + i * h));
            //}
            //int _N = 1000;
            //N = 1000;
            //h = (rangeEnd - rangeStart) / N;

            //double sum = 0;
            //for (int i = 100; i < _N; i += 2)
            //{
            //    N = i;
            //    h = (rangeEnd - rangeStart) / N;
            //    //chart1.Series[0].Points.AddXY(rangeStart + i * (rangeEnd - rangeStart) / _N, function(rangeStart + i * (rangeEnd - rangeStart) / _N));
            //    System.Diagnostics.Debug.WriteLine(i);
            //    double v = Math.Abs(9 - Integrate_MultipleSimpsons1_3());
            //    //sum += v;
            //    chart1.Series[0].Points.AddXY(i, v);
            //}
            //System.Diagnostics.Debug.WriteLine(sum);
            //System.Diagnostics.Debug.WriteLine(sum  / _N);
            //chart1.Series[0].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Point;
            //chart1.Series[1].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Area;
            ////chart1.Series[1].Color = Color.Red;
            //chart1.Series[1].Color = Color.FromArgb(255, 255, 0, 0);
            //chart1.Series[0].Color = Color.FromArgb(255, 0, 0, 255);
           Stopwatch stopwatch = new Stopwatch(); //객체 선언
                                                  //stopwatch.Start(); // 시간측정 시작
                                                  //stopwatch.Stop(); //시간측정 끝

            //System.Diagnostics.Debug.WriteLine("time : " +
            //                   stopwatch.ElapsedMilliseconds + "ms");
            //int[] ks = { 2, 3, 4, 5 };
            //int[] Ns = { 100, 1000, 10000, 100000 };
            //for(int i_k = 0; i_k < ks.Length; i_k++)
            //{
            //    for (int i_n = 0; i_n < Ns.Length; i_n++)
            //    {
            //        N = Ns[i_n];
            //        h = (rangeEnd - rangeStart) / N;
            //        double s = Integrate_RichardsonExtrapolation_MultipleSegmentTrapezoidal(ks[i_k], 4);
            //        System.Diagnostics.Debug.WriteLine("{0} {1} : {2} {3}", Ns[i_n], ks[i_k], s, Math.Abs(5.3052266004051756199648523546152094782367472194741493752946067127-s));
            //        //System.Diagnostics.Debug.WriteLine(Math.Abs(9 - s) / 0.000000045);
            //    }
            //}
            //double s = Integrate_RichardsonExtrapolation_MultipleSegmentTrapezoidal(4);
            //System.Diagnostics.Debug.WriteLine(h);
            //System.Diagnostics.Debug.WriteLine(Integrate_MultipleSimpsons1_3()
            //    - NumericalCalculus.NumericalIntegration.Integrate_MultipleSimpsons1_3(function, new NumericalCalculus.ClosedInterval(0, 3, 1000)));
            //System.Diagnostics.Debug.WriteLine(Math.Abs(9 - s) / 0.000000045);


            //stopwatch.Start();
            //System.Diagnostics.Debug.WriteLine(NumericalCalculus.NumericalIntegration.Integrate_RombergMethod_MultipleSegmentTrapezoidal(function, new NumericalCalculus.ClosedInterval(rangeStart, rangeEnd, N), 12, 12));
            //stopwatch.Stop();
            //System.Diagnostics.Debug.WriteLine(stopwatch.ElapsedMilliseconds);
            //stopwatch.Reset();
           
            //stopwatch.Start();
            //System.Diagnostics.Debug.WriteLine(NumericalCalculus.NumericalIntegration.Integrate_RombergMethod_MultipleSegmentTrapezoidal_Memo(function, new NumericalCalculus.ClosedInterval(rangeStart, rangeEnd, N), 12, 12, new Dictionary<Tuple<int, int>, double>()));
            //stopwatch.Stop();
            //System.Diagnostics.Debug.WriteLine(stopwatch.ElapsedMilliseconds);
            //stopwatch.Reset();

            //stopwatch.Start();
            //System.Diagnostics.Debug.WriteLine(NumericalCalculus.NumericalIntegration.Integrate_RombergMethod_MultipleSegmentTrapezoidal_MemoBottomUp(function, new NumericalCalculus.ClosedInterval(rangeStart, rangeEnd, N), 12, 12));
            //stopwatch.Stop();
            //System.Diagnostics.Debug.WriteLine(stopwatch.ElapsedMilliseconds);
            //stopwatch.Reset();
        }
        int N = 1000;

        double rangeStart = 0;
        double rangeEnd =3;
        double h = 1;


        void D()
        {
            chart1.Series.Add("");


            //double b_a_s = (b - a) * (b - a);

            //double A = (2*fa + 4 * fc + 2*fb);
            ////double B = (fb - fa) / (b - a) - (fb + fa) / 2 * (2 / (b - a));
            //double B = 2 * (fc - fa) / (b - a) - (3 * a + b) * A / 2;
            //double C = fa - B * a - A * a * a;

            //double c2 = (2* fa - 4* fc + 2* fb) / b_a_s;
            //double c1 = (-3*(a+3*b)* fa + 4*(a+b)* fc - 3*(3*a+b)* fb) / b_a_s;
            //double c0 = (b * (a + b)* fa - 4*a*b*fc + a*(a+b)* fb) / b_a_s;
            //double c2 = (2*x2(c) - 2*x2(a))/(b-a)  -  (3*a+ b) * c1 / 2;
            //double c3 = x2(a) - a * c2 - a * a * c1;
            // System.Diagnostics.Debug.WriteLine(c1);

            double x1 = rangeStart;
            double x3 = rangeEnd;
            double x2 = (x1+x3)/2d;

            double y1 = function2(x1);
            double y2 = function2(x2);
            double y3 = function2(x3);

            double c1 = (y3 - ((x3 * (y2 - y1) + x2 * y1 - x1 * y2) / (x2 - x1))) / (Math.Pow(x3, 2) - (x3 * (x2 + x1)) + (x1 * x2));
            double c2 = ((y2 - y1) - c1 * (Math.Pow(x2, 2) - Math.Pow(x1, 2))) / (x2 - x1);
            double c3 = y1 - (c1 * Math.Pow(x1, 2)) - (c2 * x1);
            //double denominator = (x1 - x2) * (x1 - x3) * (x2 - x3);
            //double c0 = ((y1 * (x2 - x3)) + (y2 * (x3 - x1)) + (y3 * (x1 - x2))) / denominator;
            //double c1 = ((y1 * ((x2 * x2) - (x3 * x3))) + (y2 * ((x3 * x3) - (x1 * x1))) + (y3 * ((x1 * x1) - (x2 * x2)))) / denominator;
            //double c2 = ((y1 * ((x2 * x3) - (x3 * x2))) + (y2 * ((x3 * x1) - (x1 * x3))) + (y3 * ((x1 * x2) - (x2 * x1)))) / denominator;
            //QuadraticCoefficients(x1, y1, x2, y2, x3, y3, out double c0, out double c1, out double c2);
            double f2(double f2x)
            {
                return (c1 * Math.Pow(f2x, 2)) + (c2 * f2x) + c3;
            }
            N = 1000;
            for (int i = 0; i < N; i++)
            {
                double x_i = rangeStart + (i * h);
                //chart1.Series[1].Points.AddXY(rangeStart + (i - 1) * h + 0.00000000001, function(rangeStart + i * h));
                //chart1.Series[0].Points.AddXY(x_i, (c0 * x_i * x_i) + (c1 * x_i) + c2);
                chart1.Series[1].Points.AddXY(x_i, f2(x_i));
                chart1.Series[0].Points.AddXY(x_i, function2(x_i));
            }
            chart1.Series.Add("Points");
            chart1.Series[2].Points.AddXY(x1, y1);
            chart1.Series[2].Points.AddXY(x2, y2);
            chart1.Series[2].Points.AddXY(x3, y3);
            chart1.Series[2].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Point;
            chart1.Series[2].MarkerSize = 10;
            //int _N = 10000;
            //for (int i = 0; i < _N; i++)
            //{
            //}
            chart1.Series[0].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            chart1.Series[0].BorderWidth = 5;
            chart1.Series[1].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            chart1.Series[1].BorderWidth = 5;
            //chart1.Series[1].Color = Color.Red;
            chart1.Series[1].Color = Color.FromArgb(255, 255, 0, 0);
            chart1.Series[0].Color = Color.FromArgb(255, 0, 0, 255);
        }

        double function2(double x)
        {
            return Math.Sin(x) * Math.Pow(2, x);
            return x*x;
        }

        double function(double x)
        {
            //return Math.Sin(x);
            //return x * x;
            return Math.Pow(Math.Sin(Math.Pow(2, x)), 2);
        }

        double Integrate_MultipleSegmentRectangular()
        {
            double integral = 0;
            for (int i = 0; i < N; i++)
            {
                integral += function(rangeStart + i * h) * h;
            }

            return integral;
        }
        double Integrate_MultipleSegmentTrapezoidal()
        {
            double integral = 0;
            for (int i = 0; i < N; i++)
            {
                integral += function(rangeStart + i * h) + function(rangeStart + (i + 1) * h);
            }

            return integral * h / 2;
        }

        double Integrate_MultipleSegmentTrapezoidal(double local_h)
        {
            double integral = 0;
            int local_N = (int)((rangeEnd - rangeStart) / local_h);
            for (int i = 0; i < local_N; i++)
            {
                integral += function(rangeStart + i * local_h) + function(rangeStart + (i + 1) * local_h);
            }

            return integral * local_h / 2;
        }

double Integrate_MultipleSimpsons1_3()
{
    double integral = 0;

    //for (int i = 1; i <= N / 2; i += 1)
    //{
    //    double x_i = rangeStart + (2 * i - 2) * h;
    //    double x_i_1 = rangeStart + ((2 * i) - 1) * h;
    //    integral += 2 * function(x_i);
    //    integral += 4 * function(x_i_1);
    //}
    double sum_odd = 0;
    for(int i = 1; i <= N-1; i+=2)
    {
        sum_odd += function(rangeStart + i * h);
    }
    double sum_even = 0;
    for (int i = 2; i <= N - 2; i += 2)
    {
        sum_even += function(rangeStart + i * h);
    }

    integral = 4 * sum_odd + 2 * sum_even;

    return (h / 3) *(function(rangeStart) + function(rangeEnd) + integral);
}

        double Integrate_RichardsonExtrapolation_MultipleSegmentTrapezoidal(int k, int j)
        {
            if(k == 1)
            {
                return 4/3d * Integrate_MultipleSegmentTrapezoidal(h / 2)
                    - (1 / 3d) * Integrate_MultipleSegmentTrapezoidal(h);
            }
            else
            {
                int pow_4_k = (int)Math.Pow(4, k-1);

                return (pow_4_k / (pow_4_k - 1)) * Integrate_RichardsonExtrapolation_MultipleSegmentTrapezoidal(k-1, j+1)
                    - (1 / (pow_4_k - 1)) * Integrate_RichardsonExtrapolation_MultipleSegmentTrapezoidal(k - 1, j);
            }
        }

        private void chart1_Click(object sender, EventArgs e)
        {

        }
    }
}