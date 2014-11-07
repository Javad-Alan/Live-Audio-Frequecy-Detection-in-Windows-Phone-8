using System;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Navigation;
using Microsoft.Phone.Controls;
using Microsoft.Phone.Shell;
using FrequencyTest.Resources;
using Microsoft.Xna.Framework.Audio;
using DSP;
using System.Windows.Threading;
using Microsoft.Xna.Framework;

namespace FrequencyTest
{
    public partial class MainPage : PhoneApplicationPage
    {

        private Microphone microphone = Microphone.Default;
        public byte[] buffer;

        // Constructor
        public MainPage()
        {
            InitializeComponent();

            DispatcherTimer dt = new DispatcherTimer();
            dt.Interval = TimeSpan.FromMilliseconds(50);
            dt.Tick += delegate { try { FrameworkDispatcher.Update(); } catch { } };
            dt.Start();

            microphone.BufferReady += microphone_BufferReady;
            
            
        }

        void microphone_BufferReady(object sender, EventArgs e)
        {
            
            microphone.GetData(buffer);
            

            uint sampleNumbers=Utilities.NextPowerOfTwo((uint)buffer.Length);
            double[] sampleBuffer = new double[sampleNumbers];

            int index = 0;

            for (int i = 0; i < 2048; i += 2)
            {
                sampleBuffer[index] = Convert.ToDouble(BitConverter.ToInt16((byte[])buffer, i)); index++;
            }

            double[] spectrum = FourierTransform.Spectrum(ref sampleBuffer);

            uint freqIndex=0;
            double peakFreq = FourierTransform.PeakFrequency((uint)spectrum.Length, spectrum,microphone.SampleRate, ref freqIndex);

            MessageBox.Show(peakFreq.ToString() + "hz");
            
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            if (microphone.State == MicrophoneState.Stopped)
            {
                microphone.BufferDuration = TimeSpan.FromMilliseconds(100);
                buffer = new byte[microphone.GetSampleSizeInBytes(microphone.BufferDuration)];
                microphone.Start();
                //System.Diagnostics.Debug.WriteLine("Threshold setted to:" + minimumThreshold);
            }
        }

        private void RecordStop_Click(object sender, RoutedEventArgs e)
        {
            if (microphone.State == MicrophoneState.Started)
            {
                microphone.Stop();
            }
        }

        // Sample code for building a localized ApplicationBar
        
    }
}