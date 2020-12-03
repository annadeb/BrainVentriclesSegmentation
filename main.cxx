#include <iostream> // std::cout, std::cin
#include <string> // std::string
#include <itkImage.h> // itk::Image
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageSeriesReader.h>
#include <itkImageSeriesWriter.h>
#include <itkNumericSeriesFileNames.h>
#include <itkGDCMSeriesFileNames.h> 
#include <itkGDCMImageIO.h>
#include<itkGDCMImageIO.h>
#include<itkVTKImageIO.h>
#include<itkIntensityWindowingImageFilter.h>
#include<itkThresholdImageFilter.h>
#include<itkBinaryThresholdImageFilter.h>
#include<itkOtsuMultipleThresholdsImageFilter.h>
#include<itkBinaryMorphologicalOpeningImageFilter.h>
#include<itkBinaryBallStructuringElement.h>



using PixelType = signed short;
using ImageType = itk::Image<PixelType, 2>;
using BallType = itk::BinaryBallStructuringElement<short, 2>;

using Image3DType = itk::Image<PixelType, 3>;
using ImageIOType = itk::GDCMImageIO;

using SeriesReaderType = itk::ImageSeriesReader<Image3DType>;
using Series2DWriterType = itk::ImageSeriesWriter<Image3DType, ImageType>;
using Writer3Dtype = itk::ImageFileWriter<Image3DType>;
using ReaderType = itk::ImageFileReader<ImageType>;
using WriterType = itk::ImageFileWriter<ImageType>;


using NumSeriesFileNames = itk::NumericSeriesFileNames;
int main() // glowna funkcja programu
{
try{
	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();
	ImageType::Pointer image = ImageType::New();

	SeriesReaderType::Pointer seriesReader = SeriesReaderType::New();
	Series2DWriterType::Pointer series2DWriter = Series2DWriterType::New();
	itk::GDCMSeriesFileNames::Pointer namesGen = itk::GDCMSeriesFileNames::New();
	NumSeriesFileNames::Pointer numSeriesFileNames = NumSeriesFileNames::New();
	Image3DType::Pointer image3D = Image3DType::New();
	ImageIOType::Pointer dicomIO = ImageIOType::New();
	Writer3Dtype::Pointer writer3D = Writer3Dtype::New();
	itk::GDCMImageIO::Pointer gdcmImageIO = itk::GDCMImageIO::New();


	namesGen->SetDirectory("../dane_glowa/glowa/mozg zb/Head_Neck_Standard - 153275/t2_tse_tra_5");
	itk::GDCMSeriesFileNames::SeriesUIDContainerType seriesUIDs = namesGen->GetSeriesUIDs();
	itk::GDCMSeriesFileNames::FileNamesContainerType fileNames = namesGen->GetFileNames(seriesUIDs[0]);
	for (size_t i = 0; i < seriesUIDs.size(); i++) {
		std::cout << seriesUIDs[i] << std::endl;
	}
	//Czytanie serii obrazowej 
	seriesReader->SetFileNames(namesGen->GetFileNames(seriesUIDs[0]));
	seriesReader->SetImageIO(gdcmImageIO);
	image3D= seriesReader->GetOutput();
	seriesReader->Update();


	//std::cout << image3D;
	std::cout << seriesReader->GetMetaDataDictionaryArray();

	//Zapisywanie serii obrazowej 

	numSeriesFileNames->SetSeriesFormat("..\\wyniki\\IMG\%05d.dcm");
	numSeriesFileNames->SetStartIndex(1);
	numSeriesFileNames->SetEndIndex(image3D->GetLargestPossibleRegion().GetSize()[2]);

	series2DWriter->SetFileNames(numSeriesFileNames->GetFileNames());
    series2DWriter->SetImageIO(gdcmImageIO);
	//series2DWriter->SetMetaDataDictionaryArray(seriesReader->GetMetaDataDictionaryArray());
	series2DWriter->SetInput(image3D);
	series2DWriter->Update();

	
	reader->SetFileName("../wyniki/IMG00048.dcm");
	reader->Update();
	image = reader->GetOutput();

	using FilterTypeInt = itk::IntensityWindowingImageFilter<ImageType>;
	FilterTypeInt::Pointer filter = FilterTypeInt::New();
	filter->SetInput(image);
	filter->SetOutputMaximum(100);
	filter->SetOutputMinimum(0);
	filter->SetWindowLevel(500 + 500, (500 - 500) / 2);
	writer->SetInput(filter->GetOutput());
	writer->SetFileName("../wyniki/intensywnosc_48.dcm");
	writer->Update();

	using FilterType = itk::BinaryThresholdImageFilter<ImageType, ImageType>;
	FilterType::Pointer thresholder = FilterType::New();
	thresholder->SetInput(image);
	thresholder->SetInsideValue(1);
	thresholder->SetOutsideValue(0);
	thresholder->SetLowerThreshold(600);
	thresholder->SetUpperThreshold(800);
	writer->SetInput(thresholder->GetOutput());
	writer->SetFileName("../wyniki/binaryzacja_48.dcm");
	writer->Update();

	BallType::SizeType rad;
	rad[0] = 9;
	rad[1] = 5;
	int radius = 2;
	using StructuringElementType = itk::BinaryBallStructuringElement<ImageType::PixelType, ImageType::ImageDimension>;
	StructuringElementType structuringElement;
	structuringElement.SetRadius(radius);
	structuringElement.CreateStructuringElement();

	using BinaryMorphologicalOpeningImageFilterType = itk::BinaryMorphologicalOpeningImageFilter <ImageType, ImageType, StructuringElementType>;
	BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter
		= BinaryMorphologicalOpeningImageFilterType::New();
	openingFilter->SetInput(thresholder->GetOutput());
	openingFilter->SetKernel(structuringElement);
	openingFilter->SetForegroundValue(1);
	openingFilter->Update();

	writer->SetInput(openingFilter->GetOutput());
	writer->SetFileName("../wyniki/otwarcie_48.dcm");
	writer->Update();


	//using FilterType = itk::OtsuMultipleThresholdsImageFilter<ImageType, ImageType>;
	//FilterType::Pointer thresholderOtsu = FilterType::New();
	//thresholderOtsu->SetInput(image);
	//thresholderOtsu->SetLabelOffset(0);
	//thresholderOtsu->SetNumberOfHistogramBins(100);
	//thresholderOtsu->SetNumberOfThresholds(5);

	//writer->SetInput(thresholderOtsu->GetOutput());
	//writer->SetFileName("../wyniki/intensywnosc_48.dcm");
	//writer->Update();
}
catch(itk::ExceptionObject &ex){
ex.Print(std::cout);
}
catch(std::exception &ex){
std::cout << ex.what() << std::endl;
}
catch(...){
std::cout << "Unknown error!" << std::endl;
}
std::cout << "Hit [Enter]...";
std::cin.get();
return EXIT_SUCCESS; // albo EXIT_FAILURE
}