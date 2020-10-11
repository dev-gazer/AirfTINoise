#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <QPrinter>
#include <QPrintDialog>
#include <QFontDialog>
#include <QColorDialog>
#include <QFileDialog>
#include <Eigen/Dense>

using namespace Eigen;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButtonPlot_clicked();

    void on_actionClose_triggered();

    void on_actionAbout_the_App_triggered();

    void on_actionExport_as_png_triggered();

    //Add experimental points
    void addPoint(double freq, double SPL);
    void clearData();
    void plot();

    void on_pushButtonAdd_clicked();

    void on_pushButtonClear_clicked();

    void on_actionAbout_Qt_triggered();

    void on_actionAbout_QCustomPlot_triggered();

    void on_actionExport_data_as_csv_triggered();

private:
    Ui::MainWindow *ui;
    QVector<double> xValuesQVector, yValuesQVector;
    QVector<double> qv_x, qv_y;
    VectorXd SPL, freq;

};
#endif // MAINWINDOW_H
