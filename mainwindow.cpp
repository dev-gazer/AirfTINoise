#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <Eigen/Dense>
#include <math.h>
#include <QVector>
#include <QFile>
#include <QTextStream>
#include <QFileDialog>
#include <QMessageBox>

using namespace Eigen;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    //User interface setup
    ui->setupUi(this);

    //Graph area
    ui->widgetGraph->addGraph();
    ui->widgetGraph->xAxis->setScaleType(QCPAxis::stLogarithmic);
    ui->widgetGraph->yAxis->grid()->setSubGridVisible(true);
    ui->widgetGraph->xAxis->grid()->setSubGridVisible(true);
    // give the axes some labels:
    ui->widgetGraph->xAxis->setLabel("f [Hz]");
    ui->widgetGraph->yAxis->setLabel("SPL [dBA]");
    ui->widgetGraph->xAxis->setRange(10, 20000);
    QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
    ui->widgetGraph->xAxis->setTicker(logTicker);
    ui->widgetGraph->xAxis->setNumberFormat("eb");
    ui->widgetGraph->xAxis->setNumberPrecision(0);

    //Input data
    ui->lineEditChord->setText("0.15");
    ui->lineEditChord->setValidator(new QDoubleValidator(0, 100, 2, this));
    ui->lineEditSpan->setText("0.498");
    ui->lineEditSpan->setValidator(new QDoubleValidator(0, 100, 2, this));
    ui->lineEditVelocity->setText("30");
    ui->lineEditVelocity->setValidator(new QDoubleValidator(0, 100, 2, this));
    ui->lineEditIntensity->setText("3.7");
    ui->lineEditIntensity->setValidator(new QDoubleValidator(0, 100, 2, this));
    ui->lineEditLengthScale->setText("0.0065");
    ui->lineEditLengthScale->setValidator(new QDoubleValidator(0, 100, 2, this));
    ui->lineEditAoA->setText("0");
    ui->lineEditAoA->setValidator(new QDoubleValidator(0, 100, 2, this));
    ui->lineEditDistance->setText("1.25");
    ui->lineEditDistance->setValidator(new QDoubleValidator(0, 100, 2, this));
    ui->lineEditPhi->setText("90");
    ui->lineEditPhi->setValidator(new QDoubleValidator(-180, 180, 2, this));
    ui->lineEditTheta->setText("90");
    ui->lineEditTheta->setValidator(new QDoubleValidator(-180, 180, 2, this));

    //Experimental data
    ui->widgetGraph->addGraph();
    ui->widgetGraph->graph(1)->setScatterStyle(QCPScatterStyle::ssCircle);
    ui->widgetGraph->graph(1)->setPen(QPen(Qt::red));
    ui->widgetGraph->graph(1)->setLineStyle(QCPGraph::lsNone);
    ui->doubleSpinBoxFreq->setRange(0, 100000);
    ui->doubleSpinBoxSPL->setRange(-150, 999);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButtonPlot_clicked()
{


    //Lowson
    VectorXd SPL, K, S, LFC, Aux1, Aux4, Aux5;
    double f[] = {10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000};
    double Af[] = {-70.4, -63.4, -56.7, -50.5, -44.7, -39.4, -34.6, -30.2, -26.2, -22.5, -19.1, -16.1, -13.4, -10.9, -8.6, -6.6, -4.8, -3.2, -1.9, -0.8, 0, 0.6, 1, 1.2, 1.3, 1.2, 1, 0.5, -0.1, -1.1, -2.5, -4.3, -6.6, -9.3};
    //std::vector<double> K, S, LFC, Aux1, Aux4, Aux5, SPL;
    QString CHORD=ui->lineEditChord->text();
    double c = CHORD.toDouble();
    QString SPAN=ui->lineEditSpan->text();
    double s = SPAN.toDouble();
    QString VELOCITY=ui->lineEditVelocity->text();
    double U = VELOCITY.toDouble();
    QString INTENSITY=ui->lineEditIntensity->text();
    double I = INTENSITY.toDouble()/100;
    QString LENGTH=ui->lineEditLengthScale->text();
    double L = LENGTH.toDouble();
    QString AOA=ui->lineEditAoA->text();
    double AoA = AOA.toDouble()*3.1415/180;
    QString DISTANCE=ui->lineEditDistance->text();
    double r_e = DISTANCE.toDouble();
    QString PHI=ui->lineEditPhi->text();
    double phi = PHI.toDouble()*3.1415/180;
    QString THETA=ui->lineEditTheta->text();
    double theta = THETA.toDouble()*3.1415/180;
    double rho = 1.225;
    double c0 = 340;
    double M = U/c0;
    double coef_1;
    double coef_2;
    double D_L = (pow(sin(theta),2)*(pow(sin(phi),2)))/pow((1+M*cos(theta)),4);
    double beta = pow((1-pow(M, 2)),0.5);
    double Aux = 0.5*pow(rho*c0*I*U/r_e,2)*pow(M,3)*s*L*D_L;
    if(L/c < 1){
        coef_1 = 85.95;
        coef_2 = 19/6;
    }
    else{
        coef_1 = 78.4;
        coef_2 = 7/3;
    }
    Map<RowVectorXd> freq(f, 34);
    Map<RowVectorXd> A_filter(Af, 34);
    //freq.setLinSpaced(10000, 10, 20000);
    int size_freq=freq.size();
    //int size_freq = 34;
    SPL.setZero(size_freq, 1);
    K.setZero(size_freq, 1);
    S.setZero(size_freq, 1);
    LFC.setZero(size_freq, 1);
    Aux1.setZero(size_freq, 1);
    Aux4.setZero(size_freq, 1);
    Aux5.setZero(size_freq, 1);
    for (int i=0; i<size_freq; i++){
        K(i) = 3.1415*freq(i)*c/U;
        S(i) = pow(pow((2*3.1415*K(i)/(pow(beta, 2)))+(pow((1+(2.4*K(i)/pow(beta,2))), -1)),-1), 0.5);
        LFC(i) = (1+(9*pow(AoA,2)))*10*pow(S(i),2)*M*pow(K(i),2)/(pow(beta,2));
        Aux1(i) = (10*log10(LFC(i)/(1+LFC(i)))) + coef_1;
        Aux4(i) = pow(K(i),3)/pow(1+pow(K(i),2),coef_2);
        Aux5(i) = 10*log10(Aux*Aux4(i));
        SPL(i) = (10*log10(pow(10, (Aux1(i)+Aux5(i))/10))) + A_filter(i);
    }
    // convert the Eigen objects into the std::vector form
    // .data() returns the pointer to the first memory location of the first entry of the stored object
    // https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html
    std::vector<double> xValuesStdVector(freq.data(), freq.data() + freq.rows() * freq.cols());
    std::vector<double> yValuesStdVector(SPL.data(), SPL.data() + SPL.rows() * SPL.cols());

    //convert the std::vector objects to the Qt QVector form
    QVector<double> xValuesQVector = QVector<double>::fromStdVector(xValuesStdVector);
    QVector<double> yValuesQVector = QVector<double>::fromStdVector(yValuesStdVector);




    // this is necessary for seting the axes limits
//    double x_maxValue=freq.maxCoeff();
//    double x_minValue=freq.minCoeff();

    // this is necessary for seting the axes limits
    double y_maxValue=SPL.maxCoeff();
    //double y_minValue=SPL.minCoeff();


    QCustomPlot *customPlot=ui->widgetGraph;

    // create graph and assign data to it:
    //customPlot->addGraph();
    customPlot->graph(0)->setData(xValuesQVector, yValuesQVector);
    //customPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    //customPlot->yAxis->grid()->setSubGridVisible(true);
    //customPlot->xAxis->grid()->setSubGridVisible(true);
    // give the axes some labels:
    //customPlot->xAxis->setLabel("f [Hz]");
    //customPlot->yAxis->setLabel("SPL [dBA]");
    // set axes ranges, so we see all data:
    //customPlot->xAxis->setRange(10, 20000);
    //QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
    //customPlot->xAxis->setTicker(logTicker);
    //customPlot->xAxis->setNumberFormat("eb");
    //customPlot->xAxis->setNumberPrecision(0);
    //customPlot->yAxis->setRange(y_minValue-0.1*abs(y_minValue), y_maxValue+0.1*abs(y_maxValue));
    customPlot->yAxis->setRange(0, y_maxValue+0.1*abs(y_maxValue));
    customPlot->replot();

    return;


}

void MainWindow::on_actionClose_triggered()
{
    QApplication::quit();
}

void MainWindow::on_actionAbout_the_App_triggered()
{
        QMessageBox::information(this, "Turbulent Inflow Noise Prediction", "Written by Alexandre Martuscelli Faria, from Poli-Wind and USP, based on his PhD thesis. For additional information, visit https://sites.usp.br/poli-wind/.");
}

void MainWindow::on_actionExport_as_png_triggered()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Image"), "Simulation.png", tr("Images (*.png)"));
    if (fileName != "")
        ui->widgetGraph->savePng(fileName);

}

void MainWindow::addPoint(double freq, double SPL)
{
    qv_x.append(freq);
    qv_y.append(SPL);
}

void MainWindow::clearData()
{
    qv_x.clear();
    qv_y.clear();
}

void MainWindow::plot()
{
    ui->widgetGraph->graph(1)->setData(qv_x, qv_y);
    ui->widgetGraph->replot();
    ui->widgetGraph->update();
}

void MainWindow::on_pushButtonAdd_clicked()
{
    addPoint(ui->doubleSpinBoxFreq->value(), ui->doubleSpinBoxSPL->value());
    plot();
}

void MainWindow::on_pushButtonClear_clicked()
{
    clearData();
    plot();
}

void MainWindow::on_actionAbout_Qt_triggered()
{
            QApplication::aboutQt();
}

void MainWindow::on_actionAbout_QCustomPlot_triggered()
{
    QMessageBox::information(this, "About QCustomPlot", "QCustomPlot is a Qt C++ widget for plotting and data visualization. It has no further dependencies and is well documented. This plotting library focuses on making good looking, publication quality 2D plots, graphs and charts, as well as offering high performance for realtime visualization applications. Credits to Emanuel Eichhammer. For more information, please go to: https://www.qcustomplot.com/index.php/introduction.");
}


void MainWindow::on_actionExport_data_as_csv_triggered()
{
    //Lowson
    VectorXd SPL, K, S, LFC, Aux1, Aux4, Aux5;
    double f[] = {10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000};
    double Af[] = {-70.4, -63.4, -56.7, -50.5, -44.7, -39.4, -34.6, -30.2, -26.2, -22.5, -19.1, -16.1, -13.4, -10.9, -8.6, -6.6, -4.8, -3.2, -1.9, -0.8, 0, 0.6, 1, 1.2, 1.3, 1.2, 1, 0.5, -0.1, -1.1, -2.5, -4.3, -6.6, -9.3};
    //std::vector<double> K, S, LFC, Aux1, Aux4, Aux5, SPL;
    QString CHORD=ui->lineEditChord->text();
    double c = CHORD.toDouble();
    QString SPAN=ui->lineEditSpan->text();
    double s = SPAN.toDouble();
    QString VELOCITY=ui->lineEditVelocity->text();
    double U = VELOCITY.toDouble();
    QString INTENSITY=ui->lineEditIntensity->text();
    double I = INTENSITY.toDouble()/100;
    QString LENGTH=ui->lineEditLengthScale->text();
    double L = LENGTH.toDouble();
    QString AOA=ui->lineEditAoA->text();
    double AoA = AOA.toDouble()*3.1415/180;
    QString DISTANCE=ui->lineEditDistance->text();
    double r_e = DISTANCE.toDouble();
    QString PHI=ui->lineEditPhi->text();
    double phi = PHI.toDouble()*3.1415/180;
    QString THETA=ui->lineEditTheta->text();
    double theta = THETA.toDouble()*3.1415/180;
    double rho = 1.225;
    double c0 = 340;
    double M = U/c0;
    double coef_1;
    double coef_2;
    double D_L = (pow(sin(theta),2)*(pow(sin(phi),2)))/pow((1+M*cos(theta)),4);
    double beta = pow((1-pow(M, 2)),0.5);
    double Aux = 0.5*pow(rho*c0*I*U/r_e,2)*pow(M,3)*s*L*D_L;
    if(L/c < 1){
        coef_1 = 85.95;
        coef_2 = 19/6;
    }
    else{
        coef_1 = 78.4;
        coef_2 = 7/3;
    }
    Map<RowVectorXd> freq(f, 34);
    Map<RowVectorXd> A_filter(Af, 34);
    //freq.setLinSpaced(10000, 10, 20000);
    int size_freq=freq.size();
    //int size_freq = 34;
    SPL.setZero(size_freq, 1);
    K.setZero(size_freq, 1);
    S.setZero(size_freq, 1);
    LFC.setZero(size_freq, 1);
    Aux1.setZero(size_freq, 1);
    Aux4.setZero(size_freq, 1);
    Aux5.setZero(size_freq, 1);
    for (int i=0; i<size_freq; i++){
        K(i) = 3.1415*freq(i)*c/U;
        S(i) = pow(pow((2*3.1415*K(i)/(pow(beta, 2)))+(pow((1+(2.4*K(i)/pow(beta,2))), -1)),-1), 0.5);
        LFC(i) = (1+(9*pow(AoA,2)))*10*pow(S(i),2)*M*pow(K(i),2)/(pow(beta,2));
        Aux1(i) = (10*log10(LFC(i)/(1+LFC(i)))) + coef_1;
        Aux4(i) = pow(K(i),3)/pow(1+pow(K(i),2),coef_2);
        Aux5(i) = 10*log10(Aux*Aux4(i));
        SPL(i) = (10*log10(pow(10, (Aux1(i)+Aux5(i))/10))) + A_filter(i);
    }
    QString fileName = QFileDialog::getSaveFileName(this, tr("Export as .csv"), "Simulation.csv", tr("Comma saparated values (*.csv)"));
    QFile file(fileName);
    if (file.open(QIODevice::ReadWrite)){
        QTextStream stream(&file);
        stream << "Exported simulation from AirfNoise" << Qt::endl;
        stream << "f [HZ]" << "," << "SPL [dBA]" << Qt::endl;
        for (int i = 0; i<freq.size();i++)
            stream << freq(i) << "," << SPL(i) << Qt::endl;
    }

}
