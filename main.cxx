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
#include<itkErodeObjectMorphologyImageFilter.h>
#include<itkBinaryErodeImageFilter.h>
#include<itkConnectedComponentImageFilter.h>
#include<itkGradientMagnitudeImageFilter.h>
#include <itkCropImageFilter.h>

using PixelType = signed short;
using ImageType = itk::Image<PixelType, 2>;
using BallType = itk::BinaryBallStructuringElement<short, 3>;
using InputPixelType = float;
using OutputPixelType = float;

using Image3DType = itk::Image<PixelType, 3>;
using ImageIOType = itk::GDCMImageIO;
using InputImageType = itk::Image< InputPixelType, 3 >;
using OutputImageType = itk::Image< OutputPixelType, 3 >;
using SeriesReaderType = itk::ImageSeriesReader<Image3DType>;
using Series2DWriterType = itk::ImageSeriesWriter<Image3DType, ImageType>;
using Series3DWriterType = itk::ImageSeriesWriter<Image3DType, Image3DType>;

using Writer3Dtype = itk::ImageFileWriter<Image3DType>;
using ReaderType = itk::ImageFileReader<ImageType>;
using WriterType = itk::ImageFileWriter<ImageType>;

using ConnectedComponent = itk::ConnectedComponentImageFilter<Image3DType, Image3DType, Image3DType>;
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
	Series3DWriterType::Pointer series3DWriter = Series3DWriterType::New();
	ConnectedComponent::Pointer CCImageFilter = ConnectedComponent::New();
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
	//std::cout << seriesReader->GetMetaDataDictionaryArray();

	//Zapisywanie serii obrazowej 

	//numSeriesFileNames->SetSeriesFormat("..\\wyniki\\IMG\%05d.dcm");
	//numSeriesFileNames->SetStartIndex(1);
	//numSeriesFileNames->SetEndIndex(image3D->GetLargestPossibleRegion().GetSize()[2]);

	//series2DWriter->SetFileNames(numSeriesFileNames->GetFileNames());
 //   series2DWriter->SetImageIO(gdcmImageIO);
	////series2DWriter->SetMetaDataDictionaryArray(seriesReader->GetMetaDataDictionaryArray());
	//series2DWriter->SetInput(image3D);
	//series2DWriter->Update();
	//
	series3DWriter->SetInput(image3D);
	series3DWriter->SetFileName("..\\wyniki\\img3D.vtk");
	series3DWriter->Update();
	
	//reader->SetFileName("../wyniki/IMG00048.dcm");
	//reader->Update();
	//image = reader->GetOutput();

	//przycinanie
	using CropImageFilterType = itk::CropImageFilter<Image3DType, Image3DType>;

	CropImageFilterType::Pointer cropFilter = CropImageFilterType::New();
	cropFilter->SetInput(image3D);
	// The SetBoundaryCropSize( cropSize ) method specifies the size of
	// the boundary to be cropped at both the uppper & lower ends of the
	// image eg. cropSize pixels will be removed at both upper & lower
	// extents

	Image3DType::SizeType cropSize2;
	cropSize2[0] = 90; //R
	cropSize2[1] = 87; //A
	cropSize2[2] = 42; //I

	Image3DType::SizeType cropSize3;
	cropSize3[0] = 85; //L
	cropSize3[1] = 60; //P
	cropSize3[2] = 0; //S

	cropFilter->SetUpperBoundaryCropSize(cropSize3);

	cropFilter->SetLowerBoundaryCropSize(cropSize2);

	series3DWriter->SetInput(cropFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_crop.vtk");
	series3DWriter->Update();

	//binaryzacja
	using FilterType = itk::BinaryThresholdImageFilter<Image3DType, Image3DType>;
	FilterType::Pointer thresholder = FilterType::New();
	thresholder->SetInput(cropFilter->GetOutput());
	thresholder->SetInsideValue(1);
	thresholder->SetOutsideValue(0);
	thresholder->SetLowerThreshold(600);
	
	series3DWriter->SetInput(thresholder->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_bin_ryg.vtk");
	series3DWriter->Update();


//erozja
	BallType::SizeType rad;
	rad[0] = 1;
	rad[1] = 1;
	int radius = 1;
	using StructuringElementType = itk::BinaryBallStructuringElement<ImageType::PixelType, Image3DType::ImageDimension>;
	StructuringElementType structuringElement;
	structuringElement.SetRadius(radius);
	structuringElement.CreateStructuringElement();

	using FilterErodeType = itk::BinaryErodeImageFilter<Image3DType, Image3DType, StructuringElementType>;
	FilterErodeType::Pointer erodeFilter = FilterErodeType::New();
	erodeFilter->SetInput(thresholder->GetOutput());
	erodeFilter->SetKernel(structuringElement);
	erodeFilter->SetBackgroundValue(0);
	erodeFilter->SetForegroundValue(1);
	erodeFilter->Update();

	
	series3DWriter->SetInput(erodeFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_erode_przyc.vtk");
	series3DWriter->Update();
	
	
	//otwarcie
	/*BallType::SizeType rad;
	rad[0] = 9;
	rad[1] = 5;
	int radius = 1;
	using StructuringElementType = itk::BinaryBallStructuringElement<ImageType::PixelType, Image3DType::ImageDimension>;
	StructuringElementType structuringElement;
	structuringElement.SetRadius(radius);
	structuringElement.CreateStructuringElement();

	using BinaryMorphologicalOpeningImageFilterType = itk::BinaryMorphologicalOpeningImageFilter <Image3DType, Image3DType, StructuringElementType>;
	BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter
		= BinaryMorphologicalOpeningImageFilterType::New();
	openingFilter->SetInput(thresholder->GetOutput());
	openingFilter->SetKernel(structuringElement);
	openingFilter->SetForegroundValue(1);
	openingFilter->Update();

	series3DWriter->SetInput(openingFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_op_seria.vtk");
	series3DWriter->Update();

	CCImageFilter->SetInput(image3D);
	CCImageFilter->Update();
	std::cout << CCImageFilter->GetObjectCount() << std::endl;*/
	
	
	//Magnitude

	//using FilterGradientType = itk::GradientMagnitudeImageFilter<Image3DType, Image3DType >;
	//FilterGradientType::Pointer filterG = FilterGradientType::New();
	//filterG->SetInput(image3D);
	//filterG->Update();

	//series3DWriter->SetInput(filterG->GetOutput());
	//series3DWriter->SetFileName("..\\wyniki\\img3D_mag.vtk");
	//series3DWriter->Update();


	//thresholder->SetInput(filterG->GetOutput());
	//thresholder->SetInsideValue(1);
	//thresholder->SetOutsideValue(0);
	//thresholder->SetLowerThreshold(120);
	//thresholder->SetUpperThreshold(180);

	//series3DWriter->SetInput(thresholder->GetOutput());
	//series3DWriter->SetFileName("..\\wyniki\\img3D_bin_po_mag.vtk");
	//series3DWriter->Update();




	//series3DWriter->SetInput(image3D);

	//writer->SetInput(openingFilter->GetOutput());
	//writer->SetFileName("../wyniki/otwarcie_48.dcm");
	//writer->Update();


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
