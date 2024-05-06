using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumericalCalculus
{
    public class ClosedInterval
    {
        double _Start;
        double _End;
        int _Count;
        double _Interval;

        public ClosedInterval(double Start, double End, int Count)
        {
            _Start = Start;
            _End = End;
            _Count = Count;
            ReCalculate_Interval();
        }

        public double Start
        { //캡슐화

            get { return _Start; }
            set
            {
                _Start = value;
                ReCalculate_Interval();
            }
        }
        public double End
        { //캡슐화
            get { return _End; }
            set
            {
                _End = value;
                ReCalculate_Interval();
            }
        }
        public int Count
        { //캡슐화
            get { return _Count; }
            set
            {
                _Count = value;
                ReCalculate_Interval();
            }
        }
        public double Interval
        { //캡슐화
            get { return _Interval; }
            set
            {
                _Interval = value;
                ReCalculate_Count();
            }
        }

        void ReCalculate_Interval()
        {
            _Interval = (_End - _Start) / _Count;
        }
        void ReCalculate_Count()
        {
            _Count = (int)((_End - _Start) / _Interval);
        }
    }
}
