/*
* Add the following to the .pro Qt file:
*   QT += gui
*   QT += widgets
*   QT += core network
*   LIBS += -lz
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
// #include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <complex>
#include <cmath>
#include <cstdio>

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

// Function to apply FFT on 1D sequences without external libraries
void fft(std::vector<std::complex<double>>& a) {
    int n = a.size();
    if (n <= 1) {
        return;
    }

    // Divide the vector into even and odd parts
    std::vector<std::complex<double>> even, odd;
    for (int i = 0; i < n; i += 2) {
        even.push_back(a[i]);
        odd.push_back(a[i + 1]);
    }

    // Recursively compute FFT for even and odd parts
    fft(even);
    fft(odd);

    // Combine the results
    for (int i = 0; i < n / 2; i++) {
        std::complex<double> t = std::polar(1.0, -2.0 * M_PI * i / n) * odd[i];
        a[i] = even[i] + t;
        a[i + n / 2] = even[i] - t;
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
        // |z| = sqrt(real(z)^2 + imag(z)^2)
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

// -------------- contour detection --------------

QPainterPath detectContours(const QImage& binaryImage) {
    QPainterPath contours;

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
                    contours.moveTo(x, y);
                    contours.lineTo(x + 1, y + 1);
                }
            }
        }
    }

    return contours;
}

std::vector<std::complex<double>> detectContoursV(const QImage& binaryImage) {
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

// Function to calculate the Euclidean distance between two points (complex numbers)
double distance(const Point& p1, const Point& p2) {
    // sqrt((p2.x - p1.x)^2 + (p2.y - p1.y)^2);
    return std::abs(p2 - p1);
}

// Function to connect nearby contours
void connectContours(std::vector<Point>& contours, double proximityThreshold) {
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

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);


    // ------------------- contours experiment -----------------------
    // Assuming you have a binary image (black and white)
    QImage binaryImage("bin_img.jpeg");
    QPainterPath contoursP = detectContours(binaryImage);
    // TODO: try with and without connecting the contours to see which model performs better
    std::vector<Point> contours = detectContoursV(binaryImage);
    // Set the proximity threshold (adjust as needed)
    // after a certain value (~20-30), it stops affecting the connections
    double proximityThreshold = 25.0; // this is, essentially, a hyperparameter
    // Connect nearby contours
    connectContours(contours, proximityThreshold);

    // Create a QPainterPath from the connected contours
    QPainterPath connectedPath;
    if (!contours.empty()) {
        connectedPath.moveTo(contours.front().real(), contours.front().imag());
        for (const auto& point : contours) {
            connectedPath.lineTo(point.real(), point.imag());
        }
    }

    // Create a QGraphicsScene and add the contours to it
    QGraphicsScene scene;
    scene.addPath(connectedPath);
    // scene.addPath(contoursP);
    // Create a QGraphicsView to display the scene
    QGraphicsView view(&scene);
    view.show();
    // ------------------- contours experiment -----------------------


    // qDebug() << "Downloading MNIST dataset files...";
    // DownloadDialog downloadDialog({
    //     "http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz", //train images (60 000)
    //     "http://yann.lecun.com/exdb/mnist/train-labels-idx1-ubyte.gz", //train labels
    //     "http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz",  //test images (10 000)
    //     "http://yann.lecun.com/exdb/mnist/t10k-labels-idx1-ubyte.gz"   //test labels
    // });
    // downloadDialog.show();
    // auto filePaths = downloadDialog.downloadAll();

    // const char* trainImages = "train-images-idx3-ubyte.gz";
    // const char* trainLabels = "train-labels-idx1-ubyte.gz";

    // qDebug() << "Reading MNIST train dataset...";
    // std::vector<std::vector<uint8_t>> images = readMNISTImages(trainImages);
    // std::vector<uint8_t> labels = readMNISTLabels(trainLabels);

    // for (const auto& imageBytes : images) {
    //     // Create a QImage from the vector<uint8_t>
    //     QImage image(imageBytes.data(), 28, 28, QImage::Format_Grayscale8);

    //     // Calculate contours using the custom function
    //     auto contours = detectContoursV(image); //TODO: binarize before getting contours?

    //     std::vector<double> spectrum = applyFFT(contours);

    //     // Analyze coefficients and determine how much to cut

    //     // Your analysis logic here
    //     // You may analyze the spectrum to determine which high-frequency coefficients to cut

    //     //accumulate train data somehow?
    // }
    
    // // train a model?
    
    // // test the model
    
    // // perform analysis on how much of the high frequencies we can remove 
    // // until a certain accuracy regression threshold is reached 
    // // (e.g. if the base model is 80% accurate, a 10% reduction is the cut-off point - get the difference in frequencies)

    // qDebug() << "Cleaning up MNIST dataset files...";
    // for (auto filePath : filePaths) {
    //     std::remove(filePath.c_str());
    // }
    // downloadDialog.close();
    return app.exec();
}

#include "main.moc"
