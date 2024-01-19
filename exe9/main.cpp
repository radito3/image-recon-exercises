#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QTextStream>
#include <QString>

#include <cmath>
#include <vector>

using namespace std;

const int INTENS_MIN = 0;
const int INTENS_MAX = 255;

const double EPS = 1.0e-14;

const QString file_name = "gray_naska.jpg";

void calcHisto(QImage &image, double histo[])
{
    for (int i = 0; i <= INTENS_MAX; i++)
    {
        histo[i] = 0;
    }
    for (int indx_row = 0; indx_row < image.height(); indx_row++)
    {
        quint8* ptr_row = (quint8*)(image.bits() 
                + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++)
        {
            histo[ptr_row[indx_col]]++;
        }
    }
    int numb_pix = image.height() * image.width();
    for (int i = 0; i <= INTENS_MAX; i++)
    {
        histo[i] /= numb_pix;
    }
}//calcHisto

void printHisto(double histo[])
{
    for (int i = 0; i <= INTENS_MAX; i++)
    {
        QTextStream(stdout) << i << ": " << histo[i] << Qt::endl;   
    }
}//printHisto

void thresh(QImage &image, int thr)
{
    for (int indx_row = 0; indx_row < image.height(); indx_row++)
    {
        quint8* ptr_row = (quint8*)(image.bits() 
                + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++)
        {
            ptr_row[indx_col] = (ptr_row[indx_col] < thr) ? INTENS_MIN : INTENS_MAX;
        }
    }
}//thresh

int otsu(const double histo[])
{
    // compute cumulative sums
    double p_1[INTENS_MAX + 1] = {0};
    p_1[0] = histo[0];
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        p_1[i] = p_1[i - 1] + histo[i];
    }

    // cumulative mean
    double m[INTENS_MAX + 1] = {0};
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        m[i] = m[i - 1] + i * histo[i];
    }

    // global mean
    double m_g = m[INTENS_MAX];

    // between-class
    double b_c[INTENS_MAX + 1] = {0};
    for (int i = 1; i <= INTENS_MAX; i++)
    {
        double div = (p_1[i] * (1 - p_1[i]));
        b_c[i] = 
            fabs(div < EPS) ? 0 :
            ((m_g * p_1[i] - m[i]) * (m_g * p_1[i] - m[i])) / div;
    }

    // find max
    double max = 0;
    int max_i = 0;
    for (int i = 0; i <= INTENS_MAX; i++)
    {
        if (b_c[i] > max)
        {
            max = b_c[i];
            max_i = i;
        }
    }

    QTextStream(stdout) << "Th: " << max_i << Qt::endl;

    return max_i;
}//otsu

inline double degToRad(double deg)
{
    return (M_PI / 180.0) * deg;
}//degToRad

void houghTransform(QImage &image, int numb_rho, int numb_theta)
{
    const double MIN_THETA = -90;
    const double MAX_THETA = 90;

    int diag = sqrt(image.height() * image.height() + image.width() * image.width());

    double delta_rho = 1.0 * (2 * (diag)) / numb_rho;
    double delta_theta = 180.0 / numb_theta;

    // construct accumulator space
    vector< vector<int> > accum(numb_rho, vector<int>(numb_theta));
    for (int indx_row = 0; indx_row < numb_rho; indx_row++)
    {
        for (int indx_col = 0; indx_col < numb_theta; indx_col++)
        {
            accum[indx_row][indx_col] = 0;
        }        
    }

    // for each edge (object) pixel in the image
    for (int indx_row = 0; indx_row < image.height(); indx_row++)
    {
        quint8* ptr_row = (quint8*)(image.bits() 
                + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++)
        {
            if (ptr_row[indx_col] == INTENS_MIN)
            {
                // voting procedure
                for (int indx_theta = 0; indx_theta < numb_theta; indx_theta++)
                {
                    double theta = MIN_THETA + (indx_theta * delta_theta);
                    double rho = indx_row * cos(degToRad(theta)) + indx_col * sin(degToRad(theta));
                    int indx_rho = (int)((diag) + rho) / delta_rho;
                    accum[indx_rho][indx_theta]++;
                }
            }
        }
    }

    //QTextStream(stdout) << "Diag: " << diag << endl;
}//houghTransform

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);
    QImage image;
    QLabel label;
    if (image.load(file_name))
    {
        QTextStream(stdout) << "Image loaded: " << file_name << Qt::endl;
        QTextStream(stdout) << "Format: " << image.format() << Qt::endl;

        if (image.format() == QImage::Format_Grayscale8)
        {
            // calculate histogram
            double histo[INTENS_MAX + 1];
            calcHisto(image, histo);  

            // Otsu's method to determine global threshold
            int th = otsu(histo);

            // transform image to binary
            thresh(image, th);

            // linear Hough Transform
            houghTransform(image, 180, 180);
        } 
             
        label.setPixmap(QPixmap::fromImage(image));
    }
    else
    {
        QTextStream(stdout) << "Cannot load image: " << file_name << Qt::endl;
    }

    label.show();

    return app.exec();
}
