/*
* Add the following to the .pro Qt file:
*   QT += gui
*   QT += widgets
*   QT += core network
*   LIBS += -lz
*   LIBS += -L/path/to/lib/fftw3
*   LIBS += -lfftw3
*   
* Optionally:
*   MOC_DIR = $$OUT_PWD/moc
*   OBJECTS_DIR = $$OUT_PWD/obj
*   INCLUDEPATH += $$MOC_DIR
*/
#include <QApplication>
#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>
#include <QProgressBar>
#include <QWidget>
#include <QLabel>
#include <QVBoxLayout>
#include <QProgressDialog>
#include <QFile>
#include <QEventLoop>
#include <QDebug>
// #include <QThreadPool>
#include <QImage>
#include <QPainter>
#include <QPainterPath>
#include <QGraphicsView>
#include <QGraphicsScene>

#include <initializer_list>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <vector>
#include <bitset>
#include <unordered_map>
#include <complex>
#include <cmath>
#include <cstdio>
#include <limits>
#include <cstdlib>

#include <zlib.h>
#include <fftw3.h>

class FileDownloader : public QObject {
    Q_OBJECT

public:
    FileDownloader(QUrl&& fileUrl, QProgressBar* progressBar = new QProgressBar, QObject* parent = nullptr) 
        : QObject(parent), url(fileUrl), progressBar(progressBar)
    {
        file.setFileName(url.fileName());
        if (!file.open(QIODevice::WriteOnly)) {
            qDebug() << "Error opening file for writing:" << file.errorString();
        }
        progressBar->setMinimum(0);
        progressBar->setMaximum(0);
    }

    void download() {
        if (!file.isOpen()) {
            qDebug() << "File not open. Can not download...";
            return;
        }

        auto manager = new QNetworkAccessManager(this);

        reply = manager->get(QNetworkRequest(url));

        // Connect signals/slots to handle the data as it arrives
        connect(reply, &QNetworkReply::readyRead, this, &FileDownloader::onReadyRead);
        connect(reply, &QNetworkReply::finished, this, &FileDownloader::onFinished);
        connect(reply, &QNetworkReply::downloadProgress, this, &FileDownloader::onDownloadProgress);

        // Create an event loop to wait for the request to finish
        loop.exec();
    }

    std::filesystem::path getFilePath() const {
        return file.filesystemFileName();
    }

private slots:
    void onReadyRead() {
        while (reply->bytesAvailable()) {
            // Read and write data in chunks of 32 KB
            QByteArray data = reply->read(qMin(reply->bytesAvailable(), qint64(32 * 1024)));
            qint64 bytesWritten = file.write(data);

            // Check for errors
            if (bytesWritten == -1) {
                qDebug() << "Error writing to file:" << file.errorString();
                // You might want to abort the download or handle the error accordingly
                return;
            }
        }
    }

    void onFinished() {
        // Check for errors
        if (reply->error() != QNetworkReply::NoError) {
            qDebug() << "Error downloading" << url.fileName() << ":" << reply->errorString();
        } else {
            qDebug() << "Finished downloading" << url.fileName();
        }

        // Clean up
        file.close();
        reply->deleteLater();
        loop.quit();
    }

    void onDownloadProgress(qint64 bytesReceived, qint64 bytesTotal) {
        if (!progressBar->maximum()) {
            progressBar->setMaximum(static_cast<int>(bytesTotal));
        }
        progressBar->setValue(static_cast<int>(bytesReceived));
    }

private:
    QUrl url;
    QNetworkReply* reply;
    QFile file;
    QProgressBar* progressBar;
    QEventLoop loop;
};

class DownloadDialog : public QWidget {
    Q_OBJECT

public:
    DownloadDialog(std::initializer_list<const char*> urls, QWidget* parent = nullptr) : QWidget(parent) {
        setFixedSize(400, 500);
        setWindowTitle("MNIST Downloader");

        // Create layout
        QVBoxLayout* layout = new QVBoxLayout(this);

        for (auto* url : urls) {
            // Create widgets
            auto pb = new QProgressBar(this);
            pb->setAlignment(Qt::AlignHCenter);

            QLabel* titleLabel = new QLabel(QUrl(url).fileName(), this);
            titleLabel->setAlignment(Qt::AlignHCenter);

            auto pd = new QProgressDialog;
            pd->setWindowModality(Qt::WindowModal);
            pd->setAutoClose(false);
            pd->setAutoReset(false);
            pd->setLabel(titleLabel);
            pd->setBar(pb);
            pd->setCancelButton(nullptr);

            layout->addWidget(pd);

            progressBars.emplace(url, pb);
        }
    }

    std::vector<std::filesystem::path> downloadAll() {
        std::vector<std::filesystem::path> filePaths;
        // downladed sequentially instead of concurrently because the widget
        // does not show up until after all the downloads have finished which
        // invalidates the usage of progress bars in the first place...
        for (auto [url, pb] : progressBars) {
            FileDownloader downloader(QUrl(url), pb);
            // QThreadPool::globalInstance().start([=]() { download(url, pb); });
            downloader.download();
            filePaths.push_back(downloader.getFilePath());
        }
        return filePaths;
    }

private:
    std::unordered_map<const char*, QProgressBar*> progressBars;
};

void fftWithHighFreqCut(std::vector<std::complex<double>>& input, int numHighFreqToCut = 0) {
    int size = static_cast<int>(input.size());
    fftw_complex* in = reinterpret_cast<fftw_complex*>(input.data());

    fftw_plan plan = fftw_plan_dft_1d(size, in, in, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan);

    // apply high-frequency filter (remove N highest frequencies)
    for (int i = size - numHighFreqToCut; i < size; ++i) {
        in[i][0] = 0.0; // real part
        in[i][1] = 0.0; // imaginary part
    }

    plan = fftw_plan_dft_1d(size, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);

    // return to original scaling (the inverse transform scales the elements by the size of input)
    for (auto& el : input) {
        el /= size;
    }

    fftw_destroy_plan(plan);
}

// read gzipped MNIST images
std::vector<std::vector<uint8_t>> readMNISTImages(const char* filePath) {
    gzFile file = gzopen(filePath, "rb");

    if (!file) {
        qDebug() << "Error opening file: " << filePath;
        return {};
    }

    // Read magic number and header information
    uint32_t magicNumber, numImages, numRows, numCols;
    gzread(file, &magicNumber, sizeof(magicNumber));
    gzread(file, &numImages, sizeof(numImages));
    gzread(file, &numRows, sizeof(numRows));
    gzread(file, &numCols, sizeof(numCols));

    magicNumber = __builtin_bswap32(magicNumber);
    numImages = __builtin_bswap32(numImages);
    numRows = __builtin_bswap32(numRows);
    numCols = __builtin_bswap32(numCols);

    std::vector<std::vector<uint8_t>> images(numImages, std::vector<uint8_t>(numRows * numCols));

    // Read image data
    for (uint32_t i = 0; i < numImages; ++i) {
        gzread(file, images[i].data(), numRows * numCols);
    }

    gzclose(file);

    return images;
}

// read gzipped MNIST labels
std::vector<uint8_t> readMNISTLabels(const char* filePath) {
    gzFile file = gzopen(filePath, "rb");

    if (!file) {
        qDebug() << "Error opening file:" << filePath;
        return {};
    }

    // Read magic number and header information
    uint32_t magicNumber, numLabels;
    gzread(file, &magicNumber, sizeof(magicNumber));
    gzread(file, &numLabels, sizeof(numLabels));

    magicNumber = __builtin_bswap32(magicNumber);
    numLabels = __builtin_bswap32(numLabels);

    std::vector<uint8_t> labels(numLabels);

    // Read label data
    gzread(file, labels.data(), numLabels);

    gzclose(file);

    return labels;
}

// -------------- contour detection --------------

std::vector<std::complex<double>> detectContours(const QImage& binaryImage) {
    std::vector<std::complex<double>> contours;

    for (int y = 0; y < binaryImage.height(); ++y) {
        for (int x = 0; x < binaryImage.width(); ++x) {
            // Check if the pixel is part of the object (black pixel in binary image)
            if (qRed(binaryImage.pixel(x, y)) == 0) {
                // Check neighbors to determine if it is on the boundary
                bool isBoundary = false;

                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        int nx = x + dx;
                        int ny = y + dy;

                        // Skip pixels outside the image boundary
                        if (nx >= 0 && nx < binaryImage.width() && ny >= 0 && ny < binaryImage.height()) {
                            // Check if the neighbor is outside the object (white pixel in binary image)
                            if (qRed(binaryImage.pixel(nx, ny)) == 255) {
                                isBoundary = true;
                                break;
                            }
                        }
                    }

                    if (isBoundary) break;
                }

                // If the pixel is on the boundary, add it to the contours
                if (isBoundary) {
                    contours.push_back(std::complex<double>(x, y));
                }
            }
        }
    }
    return contours;
}

// Alias for a complex number representing a point in 2D space
using Point = std::complex<double>;

// Calculate the Euclidean distance between two points (complex numbers)
double distance(const Point& p1, const Point& p2) {
    // sqrt((p2.x - p1.x)^2 + (p2.y - p1.y)^2);
    return std::abs(p2 - p1);
}

void connectContours(std::vector<Point>& contours, double proximityThreshold = 25.0) {
    std::vector<Point> connectedContours;

    for (size_t i = 0; i < contours.size(); ++i) {
        size_t j = (i + 1) % contours.size(); // Circular index to connect the last point with the first
        // Check if the endpoints of the two contours are close enough
        if (distance(contours[i], contours[j]) < proximityThreshold) {
            // Merge the two contours
            connectedContours.insert(connectedContours.end(), contours.begin() + i, contours.begin() + j + 1);
        }
    }

    contours = connectedContours;
}

// -------------- contour detection --------------

// --------------- image processing ---------------

void to_gray(QImage& image) {
    //tranformation coefficients from color to grayscale
    const double red_coef = 0.2989,
                 green_coef = 0.5870,
                 blue_coef = 0.1140;
    
    for (int row_i = 0; row_i < image.height(); row_i++) {
        QRgb* row_p = (QRgb*) image.scanLine(row_i);

        for (int col_i = 0; col_i < image.width(); col_i++) {
            auto pixel = row_p[col_i];
            int red = qRed(pixel);
            int green = qGreen(pixel);
            int blue = qBlue(pixel);

            int gray = red_coef * red + green_coef * green + blue_coef * blue;
            row_p[col_i] = qRgb(gray, gray, gray);
        }
    }
}

// very small constant to compare floating-point numbers
const double EPS = 1.0e-14;

// gray-level intensity minimum and maximum
const int INTENS_MIN = 0;
const int INTENS_MAX = 255;

// accumulate normalized histogram
void calcHisto(const QImage &image, double histo[]) {
    for (int i = 0; i <= INTENS_MAX; i++) {
        histo[i] = 0;
    }

    for (int indx_row = 0; indx_row < image.height(); indx_row++) {
        quint8* ptr_row = (quint8*)(image.bits() + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++) {
            histo[ptr_row[indx_col]]++;
        }
    }

    int numb_pix = image.height() * image.width();
    for (int i = 0; i <= INTENS_MAX; i++) {
        histo[i] /= numb_pix;
    }
}

// apply threshold on a gray-scale image to calculate binary image
void thresh(QImage &image, int thr) {
    for (int indx_row = 0; indx_row < image.height(); indx_row++) {
        quint8* ptr_row = (quint8*)(image.bits() + indx_row * image.bytesPerLine());
        for (int indx_col = 0; indx_col < image.width(); indx_col++) {
            ptr_row[indx_col] = (ptr_row[indx_col] < thr) ? INTENS_MIN : INTENS_MAX;
        }
    }
}

// Otsu's method
int otsu(const double histo[]) {
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
}

// convert a gray-scale image to binary image
void toBinImage(QImage &image) {
    if (image.format() == QImage::Format_Grayscale8) {
        double histo[INTENS_MAX + 1];
        calcHisto(image, histo);  
        int th = otsu(histo);
        thresh(image, th);
    }
}

// --------------- image processing ---------------

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    // auto contours = detectContours(binImage);
    // fftWithHighFreqCut(contours, 1);

    // QPainterPath connectedPath;
    // connectedPath.moveTo(contours.front().real(), contours.front().imag());
    // for (const auto& point : contours) {
    //     connectedPath.lineTo(point.real(), point.imag());
    // }

    // QGraphicsScene scene;
    // scene.addPath(connectedPath);

    // QGraphicsView view(&scene);
    // view.show();

    qDebug() << "Downloading MNIST dataset files...";
    DownloadDialog downloadDialog({
        "http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz", //train images (60 000)
        "http://yann.lecun.com/exdb/mnist/train-labels-idx1-ubyte.gz", //train labels
        "http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz",  //test images (10 000)
        "http://yann.lecun.com/exdb/mnist/t10k-labels-idx1-ubyte.gz"   //test labels
    });
    downloadDialog.show();
    auto filePaths = downloadDialog.downloadAll();

    const char* trainImagesFile = "train-images-idx3-ubyte.gz";
    const char* trainLabelsFile = "train-labels-idx1-ubyte.gz";
    // const char* testImagesFile = "t10k-images-idx3-ubyte.gz";
    // const char* testLabelsFile = "t10k-labels-idx1-ubyte.gz";

    qDebug() << "Reading MNIST train dataset...";
    std::vector<std::vector<uint8_t>> trainImages = readMNISTImages(trainImagesFile);
    std::vector<uint8_t> trainLabels = readMNISTLabels(trainLabelsFile);
    // std::vector<std::vector<uint8_t>> testImages = readMNISTImages(testImagesFile);
    // std::vector<uint8_t> testLabels = readMNISTLabels(testLabelsFile);

    qDebug() << "computing image contours, FFT and cutting high frequencies...";
    for (size_t i = 0; i < 10; ++i) {
        // Create a QImage from the vector<uint8_t>
        QImage image(trainImages[i].data(), 28, 28, QImage::Format_Grayscale8);
        if (!image.allGray()) {
            to_gray(image);
        }
        QImage grayImage = image.convertToFormat(QImage::Format_Grayscale8);
        toBinImage(grayImage);
        auto contours = detectContours(grayImage);
        
        // due to the low size of the images and small contours, the number of 
        // high-frequencies we can cut from the fft contour is quite low and still be
        // able to differentiate between the digits 

        fftWithHighFreqCut(contours, 1);

        std::cout << "contours for image " << i << ":";
        for (const auto& contourPoint : contours) {
            std::cout << " " << contourPoint;
        }
        std::cout << std::endl << std::endl;
    }

    qDebug() << "Cleaning up MNIST dataset files...";
    for (auto filePath : filePaths) {
        std::remove(filePath.c_str());
    }
    downloadDialog.close();
    return app.exec();
}

#include "main.moc"
