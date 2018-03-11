#!/usr/bin/env python
"""
A simple user interface for the Music 422 PAC Coder

-----------------------------------------------------------------------
© 2009 Marina Bosi & Richard E. Goldberg -- All rights reserved
-----------------------------------------------------------------------
"""

# Uses the wxPython library for the interface
import wx
import os
from pcmfile import * # to get access to PCM files
from pacfile import * # to get access to perceptually coded files


class MyFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.notebook_1 = wx.Notebook(self, -1, style=0)
        self.notebook_1_pane_2 = wx.Panel(self.notebook_1, -1)
        self.notebook_1_pane_1 = wx.Panel(self.notebook_1, -1)
        self.label_8 = wx.StaticText(self.notebook_1_pane_1, -1, "Input File")
        self.inFile = wx.TextCtrl(self.notebook_1_pane_1, -1, "")
        self.GetInputFilename = wx.Button(self.notebook_1_pane_1, -1, "...")
        self.label_9 = wx.StaticText(self.notebook_1_pane_1, -1, "Append to create output file name")
        self.outFileAppend = wx.TextCtrl(self.notebook_1_pane_1, -1, "_decoded_2.9bps")
        self.goButton = wx.Button(self.notebook_1_pane_1, -1, "Go")
        self.label_7 = wx.StaticText(self.notebook_1_pane_1, -1, "Encode Progress")
        self.encodeGauge = wx.Gauge(self.notebook_1_pane_1, -1, 100)
        self.label_6 = wx.StaticText(self.notebook_1_pane_1, -1, "Decode Progress")
        self.decodeGauge = wx.Gauge(self.notebook_1_pane_1, -1, 100)
        self.label_1 = wx.StaticText(self.notebook_1_pane_2, -1, "Number of MDCT Lines (1/2 Block)")
        self.nMDCTLines = wx.TextCtrl(self.notebook_1_pane_2, -1, "1024")
        self.label_2 = wx.StaticText(self.notebook_1_pane_2, -1, "Number of Scale Factor Bits")
        self.nScaleBits = wx.TextCtrl(self.notebook_1_pane_2, -1, "4")
        self.label_4 = wx.StaticText(self.notebook_1_pane_2, -1, "Number of Mantissa Size Bits")
        self.nMantSizeBits = wx.TextCtrl(self.notebook_1_pane_2, -1, "4")
        self.label_5 = wx.StaticText(self.notebook_1_pane_2, -1, "Target Bit Rate (bits per sample)")
        self.targetBitsPerSample = wx.TextCtrl(self.notebook_1_pane_2, -1, "2.9")

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_BUTTON, self.SetInputFile, self.GetInputFilename)
        self.Bind(wx.EVT_BUTTON, self.DoCoding, self.goButton)

    def __set_properties(self):
        self.SetTitle("PAC Coder")
        self.GetInputFilename.SetMinSize((30, 30))

    def __do_layout(self):
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_1 = wx.GridSizer(4, 2, 2, 2)
        sizer_2 = wx.BoxSizer(wx.VERTICAL)
        sizer_3 = wx.BoxSizer(wx.VERTICAL)
        sizer_5 = wx.BoxSizer(wx.VERTICAL)
        sizer_4 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_3.Add(self.label_8, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_4.Add(self.inFile, 1, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_4.Add(self.GetInputFilename, 0, wx.ALIGN_RIGHT, 0)
        sizer_3.Add(sizer_4, 0, wx.EXPAND, 0)
        sizer_5.Add(self.label_9, 0, 0, 0)
        sizer_5.Add(self.outFileAppend, 0, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_3.Add(sizer_5, 0, wx.EXPAND, 0)
        sizer_3.Add(self.goButton, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_2.Add(sizer_3, 0, wx.EXPAND, 0)
        sizer_2.Add(self.label_7, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_2.Add(self.encodeGauge, 0, wx.EXPAND, 0)
        sizer_2.Add(self.label_6, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_2.Add(self.decodeGauge, 0, wx.EXPAND, 0)
        self.notebook_1_pane_1.SetSizer(sizer_2)
        grid_sizer_1.Add(self.label_1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.nMDCTLines, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.label_2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.nScaleBits, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.label_4, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.nMantSizeBits, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.label_5, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.targetBitsPerSample, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 0)
        self.notebook_1_pane_2.SetSizer(grid_sizer_1)
        self.notebook_1.AddPage(self.notebook_1_pane_1, "Select Input File and Run")
        self.notebook_1.AddPage(self.notebook_1_pane_2, "Change Coding Parameters")
        sizer_1.Add(self.notebook_1, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()

    def SetInputFile(self,event):
        """ picks the input WAV file"""

        # file type specifier -- only allow WAV files right now
        wildcard = "WAV file (*.wav)|*.wav"

        # create the dialog
        dlg = wx.FileDialog(
            self, message="Choose an input file",
            defaultDir=os.getcwd(),
            defaultFile="",
            wildcard=wildcard,
            style=wx.FD_OPEN | wx.FD_CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # User didn't cancel, set input filename
            self.inFile.Value=dlg.GetPath()
        dlg.Destroy()


    def CheckInputs(self):
        # check that the input values are acceptable and, if so, return True
        # check inputs, give message and return False if they are not
        inFilename=self.inFile.GetValue()
        if not os.path.exists(inFilename):
            dlg = wx.MessageDialog(self, 'Input file does not exist - please check it.',
                       'Input Error!',
                       wx.OK | wx.ICON_ERROR
                       )
            dlg.ShowModal()
            dlg.Destroy()
            return False
        return True

    def DoCoding(self, event):
        """The main control of the encode/decode process"""

        # start w/ progress bars at zero percent
        self.encodeGauge.SetValue(0)
        self.decodeGauge.SetValue(0)

        # if inputs are OK, carry out coding
        if self.CheckInputs():
            self.goButton.Disable() # can't push go till done
            self.goButton.SetLabel("Coding...")

            # get info from GUI widgets
            inFilename=self.inFile.GetValue()
            outFilename=inFilename.replace(".wav",self.outFileAppend.Value+".wav")
            codeFilename=outFilename.replace(".wav",".pac")
            nMDCTLines = int(self.nMDCTLines.GetValue())
            nScaleBits = int(self.nScaleBits.GetValue())
            nMantSizeBits = int(self.nMantSizeBits.GetValue())
            targetBitsPerSample = float(self.targetBitsPerSample.GetValue())

            # encode and then decode the selected input file
            for Direction in ("Encode", "Decode"):

                # create the audio file objects
                if Direction == "Encode":
                    inFile= PCMFile(inFilename)
                    outFile = PACFile(codeFilename)
                else: # "Decode"
                    inFile = PACFile(codeFilename)
                    outFile= PCMFile(outFilename)

                # open input file
                codingParams=inFile.OpenForReading()
                nBlocks = codingParams.numSamples/nMDCTLines #roughly number of blocks to process


                # pass parameters to the output file
                if Direction == "Encode":
                    # set additional parameters that are needed for PAC file
                    # (beyond those set by the PCM file on open)
                    codingParams.nMDCTLines = nMDCTLines
                    codingParams.nScaleBits = nScaleBits
                    codingParams.nMantSizeBits = nMantSizeBits
                    codingParams.targetBitsPerSample = targetBitsPerSample
                    # tell the PCM file how large the block size is
                    codingParams.nSamplesPerBlock = codingParams.nMDCTLines

                else: # "Decode"
                    # set PCM parameters (the rest is same as set by PAC file on open)
                    codingParams.bitsPerSample = 16

                # open the output file
                outFile.OpenForWriting(codingParams) # (includes writing header)

                # Read the input file and pass its data to the output file to be written
                iBlock=0  # current block
                if Direction=="Encode": gauge=self.encodeGauge
                else: gauge=self.decodeGauge
                while True:
                    data=inFile.ReadDataBlock(codingParams)
                    if not data: break  # we hit the end of the input file
                    outFile.WriteDataBlock(data,codingParams)
                    #update progress bar
                    iBlock +=1
                    gauge.SetValue(100*iBlock/nBlocks)  # set new value
                    gauge.Refresh()     # make sure it knows to refresh
                    wx.GetApp().Yield(True)  # yields time to other events (i.e. gauge.Refresh()) waiting
                # end loop over reading/writing the blocks

                # close the files
                inFile.Close(codingParams)
                outFile.Close(codingParams)
                gauge.SetValue(100)
            # end of loop over Encode/Decode

            # we're done - give user GUI control and tell them we're done
            self.goButton.Enable() # allow access again now
            self.goButton.SetLabel("Go")
            dlg = wx.MessageDialog(self, 'File has been encoded and then decoded!',
                                   'Done Coding',
                                   wx.OK | wx.ICON_INFORMATION
                                   )
            dlg.ShowModal()
            dlg.Destroy()
        # end codingpass for OK inputs,end of function



class MyApp(wx.App):
    def OnInit(self):
        wx.InitAllImageHandlers()
        PACCoderGUI = MyFrame(None, -1, "")
        PACCoderGUI.Center(wx.BOTH)
        self.SetTopWindow(PACCoderGUI)
        PACCoderGUI.Show()
        return 1

if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()
