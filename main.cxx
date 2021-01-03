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
#include<itkConnectedThresholdImageFilter.h>
#include<itkBinaryImageToLabelMapFilter.h>
#include<itkLabelMapToLabelImageFilter.h>
#include<itkLabelStatisticsImageFilter.h>
#include<itkLabelShapeKeepNObjectsImageFilter.h>
#include<itkRescaleIntensityImageFilter.h>
#include<itkBinaryFillholeImageFilter.h>
#include<itkVotingBinaryIterativeHoleFillingImageFilter.h>
using PixelType = signed short;
using ImageType = itk::Image<PixelType, 2>;
using BallType = itk::BinaryBallStructuringElement<short, 3>;
using InputPixelType = float;
using OutputPixelType = int;

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

using ConnectedThreshold = itk::ConnectedThresholdImageFilter<Image3DType, Image3DType>;

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
	Image3DType::Pointer image3D_bin = Image3DType::New();

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

	//Image3DType::SizeType cropSize2;
	//cropSize2[0] = 90; //R
	//cropSize2[1] = 87; //A
	//cropSize2[2] = 42; //I

	//Image3DType::SizeType cropSize3;
	//cropSize3[0] = 85; //L
	//cropSize3[1] = 60; //P
	//cropSize3[2] = 0; //S

	//cropFilter->SetUpperBoundaryCropSize(cropSize3);

	//cropFilter->SetLowerBoundaryCropSize(cropSize2);

	//series3DWriter->SetInput(cropFilter->GetOutput());
	//series3DWriter->SetFileName("..\\wyniki\\img3D_crop.vtk");
	//series3DWriter->Update(); 
		//binaryzacja
	using FilterType = itk::BinaryThresholdImageFilter<Image3DType, Image3DType>;
	FilterType::Pointer thresholder = FilterType::New();
	thresholder->SetInput(image3D);
	thresholder->SetInsideValue(1);
	thresholder->SetOutsideValue(0);
	thresholder->SetLowerThreshold(300);

	series3DWriter->SetInput(thresholder->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_bin_calosc.vtk");
	series3DWriter->Update();
	image3D_bin = thresholder->GetOutput();
	//zalewanie obszaru
	using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<Image3DType, Image3DType>;

	ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
	connected->SetInput(image3D_bin);
	connected->SetBackgroundValue(1);
		connected->Update();

	std::cout << "Number of objects: " << connected->GetObjectCount() << std::endl;
	series3DWriter->SetInput(connected->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_zalane.vtk");
	series3DWriter->Update();


	//using FillholeFilterType = itk::BinaryFillholeImageFilter<Image3DType>;


	//FillholeFilterType::Pointer fillHoleFilter = FillholeFilterType::New();
	//fillHoleFilter->SetInput(image3D_bin);
	//fillHoleFilter->SetForegroundValue(1);
	//fillHoleFilter->SetFullyConnected(1);
	//fillHoleFilter->SetCoordinateTolerance(100);
	//fillHoleFilter->Update();
	//series3DWriter->SetInput(fillHoleFilter->GetOutput());
	//series3DWriter->SetFileName("..\\wyniki\\img3D_fillHole.vtk");
	//series3DWriter->Update();

	using FillholeFilterType = itk::VotingBinaryIterativeHoleFillingImageFilter<Image3DType>;
	FillholeFilterType::Pointer fillHoleFilter = FillholeFilterType::New();
	Image3DType::SizeType a;
	a[0] = 10; //R
a[1] = 10; //A
a[2] = 10; //I
	fillHoleFilter->SetInput(image3D_bin);
	fillHoleFilter->SetForegroundValue(1);
	fillHoleFilter->SetBackgroundValue(0);
	fillHoleFilter->SetRadius(a);
	fillHoleFilter->SetMaximumNumberOfIterations(1);
	fillHoleFilter->SetMajorityThreshold(2);
	//fillHoleFilter->SetFullyConnected(1);
	//fillHoleFilter->SetCoordinateTolerance(100);
	fillHoleFilter->Update();
	series3DWriter->SetInput(fillHoleFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_fillHole.vtk");
	series3DWriter->Update();
//	using LabelShapeKeepNObjectsImageFilterType = itk::LabelShapeKeepNObjectsImageFilter<OutputImageType>;
//	LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter =
//		LabelShapeKeepNObjectsImageFilterType::New();
//	labelShapeKeepNObjectsImageFilter->SetInput(connected->GetOutput());
//	labelShapeKeepNObjectsImageFilter->SetBackgroundValue(0);
//	labelShapeKeepNObjectsImageFilter->SetNumberOfObjects(1);
//	labelShapeKeepNObjectsImageFilter->SetAttribute(
//		LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
//	
//	using RescaleFilterType = itk::RescaleIntensityImageFilter<OutputImageType, Image3DType>;
//	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
//	rescaleFilter->SetOutputMinimum(0);
//	rescaleFilter->SetOutputMaximum(itk::NumericTraits<PixelType>::max());
//	rescaleFilter->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());
//
//#ifdef ENABLE_QUICKVIEW
//	QuickView viewer;
//	viewer.AddImage(image3D, true, itksys::SystemTools::GetFilenameName(argv[1]));
//
//	std::stringstream desc;
//	desc << "Largest object of " << connected->GetObjectCount() << " objects";
//	viewer.AddImage(rescaleFilter->GetOutput(), true, desc.str());
//
//	viewer.Visualize();
//#endif
//
//	//
//	using BinaryImageToLabelMapFilterType = itk::BinaryImageToLabelMapFilter<Image3DType>;
//	BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
//	binaryImageToLabelMapFilter->SetInput(image3D_bin);
//	binaryImageToLabelMapFilter->Update();
//
//	using LabelMapToLabelImageFilterType =
//		itk::LabelMapToLabelImageFilter<BinaryImageToLabelMapFilterType::OutputImageType, Image3DType>;
//	LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
//	labelMapToLabelImageFilter->SetInput(binaryImageToLabelMapFilter->GetOutput());
//	labelMapToLabelImageFilter->Update();
//
//	using LabelStatisticsImageFilterType = itk::LabelStatisticsImageFilter<Image3DType, Image3DType>;
//	LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
//	labelStatisticsImageFilter->SetLabelInput(labelMapToLabelImageFilter->GetOutput());
//	labelStatisticsImageFilter->SetInput(image3D_bin);
//	labelStatisticsImageFilter->Update();
//
//	std::cout << "Number of labels: " << labelStatisticsImageFilter->GetNumberOfLabels() << std::endl;
//	std::cout << std::endl;
//
//	using LabelPixelType = LabelStatisticsImageFilterType::LabelPixelType;
//
//	for (auto vIt = labelStatisticsImageFilter->GetValidLabelValues().begin();
//		vIt != labelStatisticsImageFilter->GetValidLabelValues().end();
//		++vIt)
//	{
//		if (labelStatisticsImageFilter->HasLabel(*vIt))
//		{
//			LabelPixelType labelValue = *vIt;
//			std::cout << "min: " << labelStatisticsImageFilter->GetMinimum(labelValue) << std::endl;
//			std::cout << "max: " << labelStatisticsImageFilter->GetMaximum(labelValue) << std::endl;
//			std::cout << "median: " << labelStatisticsImageFilter->GetMedian(labelValue) << std::endl;
//			std::cout << "mean: " << labelStatisticsImageFilter->GetMean(labelValue) << std::endl;
//			std::cout << "sigma: " << labelStatisticsImageFilter->GetSigma(labelValue) << std::endl;
//			std::cout << "variance: " << labelStatisticsImageFilter->GetVariance(labelValue) << std::endl;
//			std::cout << "sum: " << labelStatisticsImageFilter->GetSum(labelValue) << std::endl;
//			std::cout << "count: " << labelStatisticsImageFilter->GetCount(labelValue) << std::endl;
//			// std::cout << "box: " << labelStatisticsImageFilter->GetBoundingBox( labelValue ) << std::endl; // can't output
//			// a box
//			std::cout << "region: " << labelStatisticsImageFilter->GetRegion(labelValue) << std::endl;
//			std::cout << std::endl << std::endl;
//		}
//	}

//erozja
	//BallType::SizeType rad;
	//rad[0] = 1;
	//rad[1] = 1;
	//int radius = 1;
	//using StructuringElementType = itk::BinaryBallStructuringElement<ImageType::PixelType, Image3DType::ImageDimension>;
	//StructuringElementType structuringElement;
	//structuringElement.SetRadius(radius);
	//structuringElement.CreateStructuringElement();

	//using FilterErodeType = itk::BinaryErodeImageFilter<Image3DType, Image3DType, StructuringElementType>;
	//FilterErodeType::Pointer erodeFilter = FilterErodeType::New();
	//erodeFilter->SetInput(thresholder->GetOutput());
	//erodeFilter->SetKernel(structuringElement);
	//erodeFilter->SetBackgroundValue(0);
	//erodeFilter->SetForegroundValue(1);
	//erodeFilter->Update();

	//
	//series3DWriter->SetInput(erodeFilter->GetOutput());
	//series3DWriter->SetFileName("..\\wyniki\\img3D_erode_przyc.vtk");
	//series3DWriter->Update();
	
	//rozrost 
	/*ConnectedThreshold::Pointer connThres = ConnectedThreshold::New();
	connThres->SetInput(image3D);
	connThres->SetConnectivity(ConnectedThreshold::FullConnectivity);
	connThres->SetLower(0);
	connThres->SetUpper(200);
	ConnectedThreshold::IndexType seed1;
	seed1[0] = 129; seed1[1] = 136; seed1[2] = 48;
	connThres->SetSeed(seed1);
	connThres->Update();
	
	series3DWriter->SetInput(connThres->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\rozrost.vtk");
	series3DWriter->Update();*/
	///spróbuj Confidence Filter filtr medianowy 


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
