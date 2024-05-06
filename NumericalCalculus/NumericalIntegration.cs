using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericalCalculus
{
    public class NumericalIntegration
    {
        public static double Integrate_MultipleSegmentRectangular(Func<double, double> function, ClosedInterval interval)
        {
            double integral = 0;
            int N = interval.Count; double a = interval.Start; double h = interval.Interval;
            //구간 내의 모든 지점에 대해 함숫값을 미리 계산
            double[] values = FunctionCalculater.Calculate(function, interval);
            for (int i = 0; i <= N; i++)
            {
                integral += values[i];
            }

            return integral * h;
        }
        public static double Integrate_MultipleSegmentTrapezoidal(Func<double, double> function, ClosedInterval interval)
        {
            double integral = 0;
            int N = interval.Count; double a = interval.Start; double h = interval.Interval;
            //구간 내의 모든 지점에 대해 함숫값을 미리 계산
            double[] values = FunctionCalculater.Calculate(function, interval);

            for (int i = 0; i < N; i++)
            {
                integral += values[i] + values[i + 1]; //같은 지점에서의 중복된 함숫값의 연산을 제거
            }

            return integral * h / 2d;
        }
        public static double Integrate_MultipleSimpsons1_3(Func<double, double> function, ClosedInterval interval)
        {
            double integral = 0;

            int N = interval.Count;
            if (N % 2 != 0)
            {
                N--;
            }//N을 짝수로 조정
            int halfN = N / 2;
            double a = interval.Start; double h = interval.Interval;
            //구간 내의 모든 지점에 대해 함숫값을 미리 계산
            double[] values = FunctionCalculater.Calculate(function, interval);


            double sum_odd = 0; //홀수번째 지점에서의 함숫값들의 합
            for (int i = 1; i <= halfN; i++)
            {
                sum_odd += values[2 * i - 1];
            }
            double sum_even = 0; //짝수번째 지점에서의 함숫값들의 합 
            for (int i = 1; i <= halfN; i++)
            {
                sum_even += values[2 * i];
            }

            integral = 4 * sum_odd + 2 * sum_even;

            return (h / 3) * (values[0] - values[N] + integral);
        }
        public static double Integrate_RichardsonExtrapolation_MultipleSegmentTrapezoidal(Func<double, double> function, ClosedInterval interval, int k)
        {
            double pow_4_k = Math.Pow(4, k - 1);

            ClosedInterval interval1 = new ClosedInterval(interval.Start, interval.End, interval.Count);
            interval1.Interval = interval1.Interval / Math.Pow(2, k - 1);

            ClosedInterval interval2 = new ClosedInterval(interval.Start, interval.End, interval.Count);

            return (pow_4_k / (pow_4_k - 1)) * Integrate_MultipleSegmentTrapezoidal(function, interval1)
                - (1 / (pow_4_k - 1)) * Integrate_MultipleSegmentTrapezoidal(function, interval2);

        }
        public static double Integrate_RombergMethod_MultipleSegmentTrapezoidal(Func<double, double> function, ClosedInterval interval, int n, int k)
        {
            int N = interval.Count; double a = interval.Start; double h = interval.Interval;
            if (n == 1)
            {
                ClosedInterval interval1 = new ClosedInterval(interval.Start, interval.End, interval.Count);
                interval1.Interval = interval1.Interval / Math.Pow(2, k);

                ClosedInterval interval2 = new ClosedInterval(interval.Start, interval.End, interval.Count);
                interval2.Interval = interval2.Interval / Math.Pow(2, k - 1);

                return (double)(4 / 3d) * Integrate_MultipleSegmentTrapezoidal(function, interval1)
                    - (double)(1 / 3d) * Integrate_MultipleSegmentTrapezoidal(function, interval2);
            }
            else
            {
                double pow_4_n = Math.Pow(4, n);

                return (pow_4_n / (pow_4_n - 1)) * Integrate_RombergMethod_MultipleSegmentTrapezoidal(function, interval, n - 1, k)
                    - (1d / (pow_4_n - 1)) * Integrate_RombergMethod_MultipleSegmentTrapezoidal(function, interval, n - 1, k - 1);
            }
        }

        public static double Integrate_RombergMethod_MultipleSegmentTrapezoidal_Memo(Func<double, double> function, ClosedInterval interval, int n, int k, Dictionary<Tuple<int, int>, double> memoization)
        {
            // 딕셔너리에 값이 존재하면, 이미 계산된 값을 반환
            if (memoization.ContainsKey(new Tuple<int, int>(n, k)))
            {
                return memoization[new Tuple<int, int>(n, k)];
            }

            int N = interval.Count;
            double a = interval.Start;
            double h = interval.Interval;

            if (n == 1)
            {
                ClosedInterval interval1 = new ClosedInterval(interval.Start, interval.End, interval.Count);
                interval1.Interval = interval1.Interval / Math.Pow(2, k);

                ClosedInterval interval2 = new ClosedInterval(interval.Start, interval.End, interval.Count);
                interval2.Interval = interval2.Interval / Math.Pow(2, k - 1);

                double result = (4.0 / 3.0) * Integrate_MultipleSegmentTrapezoidal(function, interval1)
                                - (1.0 / 3.0) * Integrate_MultipleSegmentTrapezoidal(function, interval2);

                //값을 저장(기억)
                memoization[new Tuple<int, int>(n, k)] = result;
                return result;
            }
            else
            {
                double pow_4_n = Math.Pow(4, n);

                double result = (pow_4_n / (pow_4_n - 1)) * Integrate_RombergMethod_MultipleSegmentTrapezoidal_Memo(function, interval, n - 1, k, memoization)
                                - (1.0 / (pow_4_n - 1)) * Integrate_RombergMethod_MultipleSegmentTrapezoidal_Memo(function, interval, n - 1, k - 1, memoization);

                //값을 저장(기억)
                memoization[new Tuple<int, int>(n, k)] = result;
                return result;
            }
        }
        public static double Integrate_RombergMethod_MultipleSegmentTrapezoidal_MemoBottomUp(Func<double, double> function, ClosedInterval interval, int max_n, int max_k)
        {
            //동적 계획법을 사용할 때 일반적으로 값을 저장하는 배열은 dp라고 이름을 설정한다
            //또한, DP가 n과 k, 2개의 변수에 관한 값이므로 2차원 배열로 만들어야 한다
            //dp[n, k]의 형태로 구현
            double[,] dp = new double[max_n + 1, max_k + 1];

            // n = 1인 경우를 모두 계산
            for (int k = 1; k <= max_k; k++)
            {
                ClosedInterval interval1 = new ClosedInterval(interval.Start, interval.End, interval.Count);
                interval1.Interval = interval1.Interval / Math.Pow(2, k);
                ClosedInterval interval2 = new ClosedInterval(interval.Start, interval.End, interval.Count);
                interval2.Interval = interval2.Interval / Math.Pow(2, k - 1);

                dp[1, k] = (4.0 / 3.0) * Integrate_MultipleSegmentTrapezoidal(function, interval1)
                            - (1.0 / 3.0) * Integrate_MultipleSegmentTrapezoidal(function, interval2);
            }

            // DP 계산
            for (int n = 2; n <= max_n; n++)
            {
                for (int k = 1; k <= max_k; k++)
                {
                    double pow_4_n = Math.Pow(4, n);
                    dp[n, k] = (pow_4_n / (pow_4_n - 1)) * dp[n - 1, k]
                                - (1.0 / (pow_4_n - 1)) * dp[n - 1, k - 1];
                }
            }

            return dp[max_n, max_k];
        }
    }
}
