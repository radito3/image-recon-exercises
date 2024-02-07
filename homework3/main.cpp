/*
* Add the following to the .pro Qt file:
*   QT += gui
*   QT += widgets
*   QT += core network
*   LIBS += -lz
*   CONFIG += c++20
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
#include <QThreadPool>
#include <QImage>

#include <initializer_list>
#include <unordered_map>
// #include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#include <zlib.h>

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

    void downloadAllAsync() {
        for (auto [url, pb] : progressBars) {
            // capture by value, as these are pointers
            QThreadPool::globalInstance()->start([=]() {
                // captured structured bindings are a C++20 extension
                // an alternative is to download the files sequentially (wihout using QThreadPool)
                FileDownloader downloader(QUrl(url), pb);
                downloader.download();
            });
        }
    
        // Wait for all tasks to finish
        QThreadPool::globalInstance()->waitForDone();
    }

private:
    std::unordered_map<const char*, QProgressBar*> progressBars;
};

// Function to calculate image contours
std::vector<std::complex<double>> calculateContours(const QImage& image) {
    std::vector<std::complex<double>> contours;
    const int THRESH = 128;

    // Iterate through each pixel in the image
    for (int y = 0; y < image.height(); ++y) {
        for (int x = 0; x < image.width(); ++x) {
            // Check if the pixel is part of the contour
            if (qRed(image.pixel(x, y)) < THRESH) {
                contours.emplace_back(x, y);
            }
        }
    }

    return contours;
}

// Function to apply FFT on 1D sequences without external libraries
void fft(std::vector<std::complex<double>> &a) {
    int n = a.size();
    if (n <= 1) {
        return;
    }

    // Divide
    std::vector<std::complex<double>> a0(a.begin(), a.begin() + n / 2);
    std::vector<std::complex<double>> a1(a.begin() + n / 2, a.end());

    // Calculate halves
    fft(a0);
    fft(a1);

    // Combine
    for (int i = 0; i < n / 2; ++i) {
        std::complex<double> t = std::polar(1.0, -2.0 * M_PI * i / n) * a1[i];
        a[i] = a0[i] + t;
        a[i + n / 2] = a0[i] - t;
    }
}

// Function to apply FFT on 1D sequences and return magnitude spectrum
std::vector<double> applyFFT(const std::vector<std::complex<double>> &contours) {
    int N = contours.size();
    std::vector<std::complex<double>> input(contours);

    fft(input);

    // Extract magnitude spectrum
    std::vector<double> spectrum(N);
    for (int i = 0; i < N; i++) {
        spectrum[i] = std::abs(input[i]);
    }

    return spectrum;
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

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    DownloadDialog downloadDialog({
        "http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz", //train imgs
        "http://yann.lecun.com/exdb/mnist/train-labels-idx1-ubyte.gz", //train labels
        "http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz",  //test img
        "http://yann.lecun.com/exdb/mnist/t10k-labels-idx1-ubyte.gz"   //test labels
    });
    downloadDialog.show();
    downloadDialog.downloadAllAsync();

    // Load MNIST dataset image (assuming a single digit image for simplicity)
    const char* trainImages = "train-images-idx3-ubyte.gz";
    const char* trainLabels = "train-labels-idx1-ubyte.gz";

    // Read gzipped MNIST images and labels
    std::vector<std::vector<uint8_t>> images = readMNISTImages(trainImages);
    std::vector<uint8_t> labels = readMNISTLabels(trainLabels);

    for (const auto& imageBytes : images) {
        // Create a QImage from the vector<uint8_t>
        QImage image(imageBytes.data(), 28, 28, QImage::Format_Grayscale8);

        // Calculate contours using the custom function
        auto contours = calculateContours(image); //TODO: binarize before getting contours? (moore.cpp)

        std::vector<double> spectrum = applyFFT(contours);

        // Analyze coefficients and determine how much to cut

        // Your analysis logic here
        // You may analyze the spectrum to determine which high-frequency coefficients to cut

        //accumulate train data somehow?
    }
    
    // train a model?
    
    // test the model

    return app.exec();
}

#include "main.moc"
