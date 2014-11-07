Live-Audio-Frequecy-Detection-in-Windows-Phone-8
================================================

Code Sample -  Windows Phone 8 Live Audio Frequency Detection

This is a simple code sample that shows how to detect audio frequency
from microphone buffer fetched from a live recording in Windows Phone 8.

###Step 1:
The first step would be getting microphone buffer data. To do that you can
go ahead and use Microphone Class and use it. 

```
private Microphone microphone = Microphone.Default;
public MainPage()
{
    InitializeComponent();

    DispatcherTimer dt = new DispatcherTimer();
    dt.Interval = TimeSpan.FromMilliseconds(50);
    dt.Tick += delegate { try { FrameworkDispatcher.Update(); } catch { } };
    dt.Start();

    microphone.BufferReady += microphone_BufferReady;
               
}
```

The thing that catches eye first is the DispatcherTimer usage, and definitely 
using the FrameworkDispatcher.Update() method that comes off the xna family.
The reason is pretty understandable when you look off to the microphone event 
we are firing here that is named BufferReady. That gets invoked when buffer is 
ready and we need to clear off the buffer and update the whole framework so.

###Step 2:
Getting a Fast Fourier Transformation going would be the next thing to do as it 
would lead you to land of frequency domain from time domain. FFT is essentially
discrete fourier transformation an as the name sounds it needs a discrete signal.
As sound is continuous we can sample our microphone feed into a discreet signal.
In this case we sampled it on 16000 hz based on the microphone.SampleRate property.

Then the nyquist or folding frequency comes in play. For a signal sampled with the 
sampling rate fs, due to the limitation imposed by the Shannon Theorem, components 
of the spectrum can only be extracted for frequencies between -fs/2 <= f <= fs/2, 
or between the negative and positive Nyquist frequency.  

To read more you guys definitely can check these awesome wikis:

(http://developer.nokia.com/community/wiki/Sound_pattern_matching_using_Fast_Fourier_Transform_in_Windows_Phone)
(http://developer.nokia.com/community/wiki/How_to_access_and_manage_the_Microphone_raw_data_in_WP)

Now here, with the new DSP class provided the FFT and frequency detection became a
breeze.

```
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

            Debug.Writeline(peakFreq.ToString() + "hz");
            
        }
```

Now, as per Nyquist sample theorem goes our example is suitable for only detecting 8Khz frequency. Each sample is 2 bytes
long and 32 KB in size. 

Will try to find a way to do this on Windows Phone 8.1


