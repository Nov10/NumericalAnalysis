using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericalCalculus {
    public class FunctionCalculater {
        public static double[] Calculate(Func<double, double> function, ClosedInterval interval)
        {
            int N = interval.Count;
            double a = interval.Start;
            double h = interval.Interval;

            double[] result = new double[N+1];
            for (int i = 0; i <= N; i++){
                result[i] = function(a + i * h);
            }

            return result;
        }
    }
}
