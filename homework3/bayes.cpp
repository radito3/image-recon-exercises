#include <QApplication>
#include <QImage>
#include <QLabel>
#include <QTextStream>
#include <QString>

#include <utility>
#include <vector>
#include <cmath>

using namespace std;

// training data type
// first dimension is feature id, second is sequence of <feature value, class id>
typedef vector< vector< pair<int, int> > > TrainData;

// very small constant to compare floating-point numbers
const double EPS = 1.0e-14;

// gray-level intensity minimum and maximum
const int INTENS_MIN = 0;
const int INTENS_MAX = 255;

// number of classes: background and object
const int CLS_NUMB = 2;

// number of features: red, green, blue
const int FTR_NUMB = 3;

// input images file names
const QString file_name_trn = "../manuscript1.jpg"; // training dataset
const QString file_name_tst = "../manuscript2.jpg"; // testing dataset

// convert color (r, g, b) pixels to equal gray-scale pixels 
void toGray(QImage &image)
{
    const double RED_COEF = 0.2989;
    const double GRN_COEF = 0.5870;
    const double BLU_COEF = 0.1140; 
    for (int indx_row = 0; indx_row < image.height(); indx_row++)
    {
        QRgb* pnt_row = (QRgb*)image.scanLine(indx_row);
        for (int indx_col = 0; indx_col < image.width(); indx_col++)
        {
            int red = qRed(pnt_row[indx_col]);
            int grn = qGreen(pnt_row[indx_col]);
            int blu = qBlue(pnt_row[indx_col]);
            int grey = RED_COEF * red + GRN_COEF * grn + BLU_COEF * blu;
            pnt_row[indx_col] = qRgb(grey, grey, grey);
                               
        }
    }
}//toGray

// convert a color image to gray-scale image
void toGrayImage(QImage &image)
{
    // covert pixels
    toGray(image);

    if (image.allGray())
    {
        // convert format of the input image to gray-scale
        image = image.convertToFormat(QImage::Format_Grayscale8);
    }
}//toGrayImage

// accumulate normalized histogram
void calcHisto(const QImage &image, double histo[])
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

// apply threshold on a gray-scale image to calculate binary image
void thresh(QImage &image, int thr)
{
    for (int indx_row = 0; indx_row < image.height(); indx_row++)
    {
        quint8* ptr_row = (quint8*)(image.bits() 
                + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++)
        {
            ptr_row[indx_col] = 
                (ptr_row[indx_col] < thr) ? INTENS_MIN : INTENS_MAX;
        }
    }
}//thresh

// Otsu's method
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

    return max_i;
}//otsu

// convert a gray-scale image to binary image
void toBinImage(QImage &image)
{
    if (image.format() == QImage::Format_Grayscale8)
    {
        double histo[INTENS_MAX + 1];
        calcHisto(image, histo);  
        int th = otsu(histo);
        thresh(image, th);
    }
}//toBinImage

// training data extraction
void trainDataExtract(
        const QImage &image_org,    // original color image 
        const QImage &image_prc,    // processed binary image with pixels labeled as background / object
        TrainData &data             // extracted training data
        )
{
    // resize to the number of features: red, green, blue
    data.resize(FTR_NUMB);

    for (int indx_row = 0; indx_row < image_org.height(); indx_row++)
    {
        QRgb* pnt_row_org = (QRgb*)image_org.scanLine(indx_row);
        quint8* ptr_row_prc = (quint8*)(image_prc.bits() 
                + indx_row * image_prc.bytesPerLine());
        for (int indx_col = 0; indx_col < image_org.width(); indx_col++)
        {
            // class id from binary image: 0 - background, 1 - object  
            int class_id = (ptr_row_prc[indx_col] == INTENS_MIN) ? 1 : 0;

            // red feature
            int red = qRed(pnt_row_org[indx_col]);
            data[0].push_back(pair<int, int>(red, class_id));

            // green feature
            int grn = qGreen(pnt_row_org[indx_col]);
            data[1].push_back(pair<int, int>(grn, class_id));

            // blue feature
            int blu = qBlue(pnt_row_org[indx_col]);
            data[2].push_back(pair<int, int>(blu, class_id));                                          
        }
    }
}//trainDataExtract

// train the model, accumulate class probability and pdf
void train(
        const TrainData &data,
        vector<double> &class_pr,
        vector< vector<double> > &pdf0,
        vector< vector<double> > &pdf1
        )
{
    // initialize probability vectors
    class_pr[0] = class_pr[1] = 0;
    
    for (int indx_ftr = 0; indx_ftr < FTR_NUMB; indx_ftr++)
    {
        pdf0[indx_ftr].resize(INTENS_MAX + 1); // resize in [0, 255]
        pdf1[indx_ftr].resize(INTENS_MAX + 1); // resize in [0, 255]
        for (int indx_lvls = 0; indx_lvls <= INTENS_MAX; indx_lvls++)
        {
            pdf0[indx_ftr][indx_lvls] = 0;
            pdf1[indx_ftr][indx_lvls] = 0;
        }
    }

    // accumulate class probabilities
    // its enough to accumulate through one of the features because each pixel is (r, g, b)   
    for (unsigned long indx_rec = 0; indx_rec < data[0].size(); indx_rec++)
    {
        class_pr[data[0][indx_rec].second]++;
    }       

    class_pr[0] /= data[0].size();
    class_pr[1] /= data[0].size();

    // accumulate probability density functions
    for (int indx_ftr = 0; indx_ftr < FTR_NUMB; indx_ftr++)
    {
        for (unsigned long indx_rec = 0; indx_rec < data[indx_ftr].size(); indx_rec++) 
        {
            data[indx_ftr][indx_rec].second ?
            pdf1[indx_ftr][data[indx_ftr][indx_rec].first]++ :
            pdf0[indx_ftr][data[indx_ftr][indx_rec].first]++;
        }
        for (int indx_lvls = 0; indx_lvls <= INTENS_MAX; indx_lvls++)
        {
                pdf0[indx_ftr][indx_lvls] /= data[indx_ftr].size();
                pdf1[indx_ftr][indx_lvls] /= data[indx_ftr].size();
                //QTextStream(stdout) << pdf0[indx_ftr][indx_lvls] << " " << pdf1[indx_ftr][indx_lvls] <<  endl;
        }
    }
    //QTextStream(stdout) << class_pr[0] << " " << class_pr[1] << endl;
}//train

void test(
        QImage &image_tst,
        const vector<double> &class_pr,
        const vector< vector<double> > &pdf0,
        const vector< vector<double> > &pdf1
        )
{
    for (int indx_row = 0; indx_row < image_tst.height(); indx_row++)
    {
        QRgb* pnt_row = (QRgb*)image_tst.scanLine(indx_row);
        for (int indx_col = 0; indx_col < image_tst.width(); indx_col++)
        {
            int red = qRed(pnt_row[indx_col]);
            int grn = qGreen(pnt_row[indx_col]);
            int blu = qBlue(pnt_row[indx_col]);

            double vote_0 = class_pr[0] * pdf0[0][red] * pdf0[1][grn] * pdf0[2][blu];
            double vote_1 = class_pr[1] * pdf1[0][red] * pdf1[1][grn] * pdf1[2][blu];

            int value = vote_0 > vote_1 ? INTENS_MAX : INTENS_MIN;
            
            pnt_row[indx_col] = qRgb(value, value, value);
                               
        }
    }
}//test

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);
    QImage image_org;       // original color image for training
    QImage image_prc;       // image to be converted to binary for training
    QImage image_tst;       // test image
    QLabel label;
    if (image_org.load(file_name_trn))
    {
        QTextStream(stdout) << "Image loaded: " << file_name_trn << endl;
        QTextStream(stdout) << "Format: " << image_org.format() << endl;

        // copy image training to be processed
        image_prc = image_org;

        if (image_prc.format() == QImage::Format_RGB32)
        {
            // to gray-scale image
            toGrayImage(image_prc);

            // convert to binary image
            toBinImage(image_prc);

            // extract training data
            TrainData data;             // data structure to store training data
            trainDataExtract(image_org, image_prc, data);

            // probabilities vectors that comprise the model
            vector<double> class_pr(CLS_NUMB);          // class probabilities for 0 - background, 1 - object
            vector< vector<double> > pdf0(FTR_NUMB);    // probability density function for background
            vector< vector<double> > pdf1(FTR_NUMB);    // probability density function for object

            // train the model
            train(data, class_pr, pdf0, pdf1);

            image_tst.load(file_name_tst);

            // test the model
            test(image_tst, class_pr, pdf0, pdf1);
        } 
    
        label.setPixmap(QPixmap::fromImage(image_tst));
        image_tst.save("out.png");
    }
    else
    {
        QTextStream(stdout) << "Cannot load all images. " << endl;
    }

    label.show();

    return app.exec();
}//main
