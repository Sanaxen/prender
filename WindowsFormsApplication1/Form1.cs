using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace WindowsFormsApplication1
{
    public partial class Form1 : Form
    {
        public string targetDir = "";
        public int convPng = 0;
        public int autoclose = 0;

        int retry;
        int count2;
        int count;
        int step;
        public Form1()
        {
            count = 0;
            count2 = 1;
            retry = 0;
            step = 0;
            InitializeComponent();
        }

        private bool draw()
        {
            if (step == 0) return false;
            string file = targetDir + "\\" + string.Format("image\\output_W{0,0:D6}.bmp", count);
            string file2 = targetDir + "\\" + string.Format("image\\output_W{0,0:D6}.bmp", count + step);

            string f = targetDir + "\\" + string.Format("output\\Closed_{0,0:D6}", count);
            int endFlg = 0;
            if (System.IO.File.Exists(f) && checkBox3.Checked )
            {
                endFlg = 1;
                //Close();
            }



            if (System.IO.File.Exists(file) || endFlg == 1)
            {
                if (endFlg == 0 && !System.IO.File.Exists(file2))
                {
                    retry++;
                    if (retry * timer1.Interval / 1000 < 60 * 3)
                    {
                        return false;
                    }
                }

                retry = 0;
                Bitmap bmp = new Bitmap(file);
                if ( true )
                {
                    this.pictureBox1.Image = bmp;

                    pictureBox1.Refresh();


                    Bitmap saveImg = null;
                    Graphics gg = null;


                    if (gg != null)
                    {
                        string savefile = targetDir + "\\" + string.Format("image\\output_W{0,0:D6}.png", count);

                        saveImg.Save(savefile, System.Drawing.Imaging.ImageFormat.Png);
                        saveImg.Dispose();
                        gg.Dispose();
                    }

                    if (count / step >= this.hScrollBar1.Maximum)
                    {
                        this.hScrollBar1.Maximum += 1000;
                    }
                    this.hScrollBar1.Value = count / step;

                    if (endFlg == 1)
                    {
                        Close();
                    }
                    return true;
                }
            }
            return false;
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            if (step == 0) return;

            if (draw())
            {
                count += step;
                if (count / step >= this.hScrollBar1.Maximum)
                {
                    this.hScrollBar1.Maximum += 1000;
                }
                this.hScrollBar1.Value = count / step;
            }
        }

        private void timer2_Tick(object sender, EventArgs e)
        {
            if (step != 0) return;

            string file = targetDir + "\\" + string.Format("image\\output_W{0,0:D6}.bmp", 0);


            if (System.IO.File.Exists(file))
            {
                string file2 = targetDir + "\\" + string.Format("image\\output_W{0,0:D6}.bmp", count2);
                if (System.IO.File.Exists(file2))
                {
                    step = count2;
                }
                else
                {
                    count2++;
                    if (count2 > 1000) count2 = 1;
                }
            }
        }

        private void hScrollBar1_ValueChanged(object sender, EventArgs e)
        {
            count = this.hScrollBar1.Value * step;
        }

        private void button2_Click(object sender, EventArgs e)
        {
            timer1.Enabled = false;
        }

        private void button3_Click(object sender, EventArgs e)
        {
            timer1.Enabled = true;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            count = 0;
            this.hScrollBar1.Value = 0;
        }

        private void hScrollBar1_Scroll(object sender, ScrollEventArgs e)
        {
            draw();
        }

        private void Form1_ResizeBegin(object sender, EventArgs e)
        {
        }

        private void Form1_Resize(object sender, EventArgs e)
        {
        }

        private void Form1_ResizeEnd(object sender, EventArgs e)
        {
        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {

        }

        private void checkBox2_CheckedChanged(object sender, EventArgs e)
        {

        }
    }
}
