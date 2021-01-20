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
#include<itkBinaryBallStructuringElement.h>
#include<itkErodeObjectMorphologyImageFilter.h>
#include<itkBinaryErodeImageFilter.h>
#include<itkBinaryDilateImageFilter.h>
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
#include<itkMultiplyImageFilter.h>
#include<itkMedianImageFilter.h>
#include<itkMinimumMaximumImageCalculator.h>
#include<itkBinaryMorphologicalOpeningImageFilter.h>
#include<itkConfidenceConnectedImageFilter.h>
//#include<itkCurvatureFlowImageFilter.h>
#include<itkNeighborhoodConnectedImageFilter.h>
//#include <itkAbsoluteValueDifferenceImageFilter.h>
#include <itkInvertIntensityImageFilter.h>
#include<itkRelabelComponentImageFilter.h>
#include<itkFlipImageFilter.h>
#include<itkAddImageFilter.h>
#include<itkJoinImageFilter.h>
#include "itkLabelGeometryImageFilter.h"
#include<itkImageRegistrationMethod.h>





using PixelType = signed short;
using ImageType = itk::Image<PixelType, 2>;
using BallType = itk::BinaryBallStructuringElement<short, 3>;
//using InputPixelType = float;
//using OutputPixelType = int;

using Image3DType = itk::Image<PixelType, 3>;
using ImageIOType = itk::GDCMImageIO;
//using InputImageType = itk::Image< InputPixelType, 3 >;
//using OutputImageType = itk::Image< OutputPixelType, 3 >;
using SeriesReaderType = itk::ImageSeriesReader<Image3DType>;
using Series2DWriterType = itk::ImageSeriesWriter<Image3DType, ImageType>;
using Series3DWriterType = itk::ImageSeriesWriter<Image3DType, Image3DType>;

using Writer3Dtype = itk::ImageFileWriter<Image3DType>;
using ReaderType = itk::ImageFileReader<ImageType>; 
using Reader3DType = itk::ImageFileReader<Image3DType>;

using WriterType = itk::ImageFileWriter<ImageType>;

using ConnectedComponent = itk::ConnectedComponentImageFilter<Image3DType, Image3DType, Image3DType>;
using NumSeriesFileNames = itk::NumericSeriesFileNames;

using ConnectedThreshold = itk::ConnectedThresholdImageFilter<Image3DType, Image3DType>;

int main() // glowna funkcja programu
{
	try {
		ReaderType::Pointer reader = ReaderType::New();
		Reader3DType::Pointer reader3D = Reader3DType::New();
		Reader3DType::Pointer reader3Db = Reader3DType::New();
		Reader3DType::Pointer reader3Dc = Reader3DType::New();
		Reader3DType::Pointer reader3Dd = Reader3DType::New();
		WriterType::Pointer writer = WriterType::New();
		ImageType::Pointer image = ImageType::New();

		SeriesReaderType::Pointer seriesReader = SeriesReaderType::New();
		Series2DWriterType::Pointer series2DWriter = Series2DWriterType::New();
		itk::GDCMSeriesFileNames::Pointer namesGen = itk::GDCMSeriesFileNames::New();
		NumSeriesFileNames::Pointer numSeriesFileNames = NumSeriesFileNames::New();
		Image3DType::Pointer image3D = Image3DType::New();
		Image3DType::Pointer image3D_bin = Image3DType::New();
		Image3DType::Pointer image3D_mnozenie = Image3DType::New();

		ImageIOType::Pointer dicomIO = ImageIOType::New();
		Writer3Dtype::Pointer writer3D = Writer3Dtype::New();
		Series3DWriterType::Pointer series3DWriter = Series3DWriterType::New();
		ConnectedComponent::Pointer CCImageFilter = ConnectedComponent::New();
		itk::GDCMImageIO::Pointer gdcmImageIO = itk::GDCMImageIO::New();


		//namesGen->SetDirectory("../dane_glowa/glowa/mozg md/Head_Neck_Standard - 1/t2_tse_tra_6");
		//namesGen->SetDirectory("../dane_glowa/glowa/mozg sd/Head_Neck_Standard - 1/t2_tse_tra_6");
		
		namesGen->SetDirectory("../dane_glowa/glowa/mozg zb/Head_Neck_Standard - 153275/t2_tse_tra_5");
		itk::GDCMSeriesFileNames::SeriesUIDContainerType seriesUIDs = namesGen->GetSeriesUIDs();
		itk::GDCMSeriesFileNames::FileNamesContainerType fileNames = namesGen->GetFileNames(seriesUIDs[0]);
		for (size_t i = 0; i < seriesUIDs.size(); i++) {
			std::cout << seriesUIDs[i] << std::endl;
		}
		//Czytanie serii obrazowej 
		seriesReader->SetFileNames(namesGen->GetFileNames(seriesUIDs[0]));
		seriesReader->SetImageIO(gdcmImageIO);
		image3D = seriesReader->GetOutput();
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
		series3DWriter->SetFileName("..\\wyniki\\1_img3D.vtk");
		series3DWriter->Update();


		//binaryzacja
		using FilterType = itk::BinaryThresholdImageFilter<Image3DType, Image3DType>;
		FilterType::Pointer thresholder = FilterType::New();
		thresholder->SetInput(image3D);
		thresholder->SetInsideValue(1);
		thresholder->SetOutsideValue(0);
		thresholder->SetLowerThreshold(300);

		series3DWriter->SetInput(thresholder->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\2_img3D_bin_calosc.vtk");
		series3DWriter->Update();
		image3D_bin = thresholder->GetOutput();


		//dylatacja
		BallType::SizeType rad;
		rad[0] = 8;
		rad[1] = 8;
		rad[2] = 8;
		using StructuringElementType = itk::BinaryBallStructuringElement<ImageType::PixelType, Image3DType::ImageDimension>;
		StructuringElementType structuringElement;
		structuringElement.SetRadius(rad);
		structuringElement.CreateStructuringElement();

		using FilterDilateType = itk::BinaryDilateImageFilter<Image3DType, Image3DType, StructuringElementType>;
		FilterDilateType::Pointer dilateFilter = FilterDilateType::New();
		dilateFilter->SetInput(image3D_bin);
		dilateFilter->SetKernel(structuringElement);
		dilateFilter->SetBackgroundValue(0);
		dilateFilter->SetForegroundValue(1);
		dilateFilter->Update();

		series3DWriter->SetInput(dilateFilter->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\3_img3D_dilate.vtk");
		series3DWriter->Update();

		//erozja 
		rad[0] = 38;
		rad[1] = 38;
		rad[2] = 38;
		structuringElement.SetRadius(rad);
		structuringElement.CreateStructuringElement();

		using FilterErodeType = itk::BinaryErodeImageFilter<Image3DType, Image3DType, StructuringElementType>;
		FilterErodeType::Pointer erodeFilter = FilterErodeType::New();
		erodeFilter->SetInput(dilateFilter->GetOutput());
		erodeFilter->SetKernel(structuringElement);
		erodeFilter->SetBackgroundValue(0);
		erodeFilter->SetForegroundValue(1);
		erodeFilter->Update();

		series3DWriter->SetInput(erodeFilter->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\4_img3D_erozja.vtk");
		series3DWriter->Update();

		//mno¿enie
		//image3D* erodeFilter->GetOutput();
		using MultiplyType = itk::MultiplyImageFilter< Image3DType, Image3DType, Image3DType>;
		MultiplyType::Pointer multiply = MultiplyType::New();
		multiply->SetInput1(image3D);
		multiply->SetInput2(erodeFilter->GetOutput());
		image3D_mnozenie = multiply->GetOutput();;
		multiply->Update();
		series3DWriter->SetInput(image3D_mnozenie);
		series3DWriter->SetFileName("..\\wyniki\\5_img3D_mnozenie.vtk");
		series3DWriter->Update();

		//filter medianowy 
		//using MedianFilterType = itk::MedianImageFilter<Image3DType, Image3DType>;
		//MedianFilterType::Pointer medianFilter = MedianFilterType::New();
		//MedianFilterType::InputSizeType radius;
		//radius.Fill(1);
		//medianFilter->SetRadius(radius);
		//medianFilter->SetInput(multiply->GetOutput());	//using FillholeFilterType = itk::BinaryFillholeImageFilter<Image3DType>;
		//series3DWriter->SetInput(medianFilter->GetOutput());
		//series3DWriter->SetFileName("..\\wyniki\\img3D_filtr_medianowy.vtk");
		//series3DWriter->Update();

		//binaryzacja -> poszukiwanie max
		thresholder->SetInput(multiply->GetOutput());
		thresholder->SetInsideValue(1);
		thresholder->SetOutsideValue(0);
		thresholder->SetLowerThreshold(600);
		thresholder->SetUpperThreshold(800);
		series3DWriter->SetInput(thresholder->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\6_img3D_bin_komory.vtk");
		series3DWriter->Update();

		//otwarcie
		using OpeningFilterType = itk::BinaryMorphologicalOpeningImageFilter<Image3DType, Image3DType, StructuringElementType>;
		int radiusOpen = 1;
		structuringElement.SetRadius(radiusOpen);
		structuringElement.CreateStructuringElement();

		OpeningFilterType::Pointer openFilter = OpeningFilterType::New();
		openFilter->SetInput(thresholder->GetOutput());
		openFilter->SetKernel(structuringElement);
		openFilter->SetBackgroundValue(0);
		openFilter->SetForegroundValue(1);
		openFilter->Update();

		series3DWriter->SetInput(openFilter->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\7_img3D_otwarcie.vtk");
		series3DWriter->Update();

		//pierwszy niezerowy punkt obrazu od czubka g³owy w dó³
		//TODO

		//rozrost
	/*	using CurvatureFlowImageFilterType =itk::CurvatureFlowImageFilter< Image3DType, Image3DType >;
		CurvatureFlowImageFilterType::Pointer smoothing =CurvatureFlowImageFilterType::New();
		smoothing->SetInput(image3D);
		smoothing->SetNumberOfIterations(5);
		smoothing->SetTimeStep(0.125);*/

		using ConnCompFilterType = itk::ConnectedComponentImageFilter<Image3DType, Image3DType>;

		ConnCompFilterType::Pointer connComp1 = ConnCompFilterType::New();
		connComp1->SetInput(openFilter->GetOutput());
		connComp1->SetBackgroundValue(0);
		connComp1->Update();
		series3DWriter->SetInput(connComp1->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\8_img3D_etykiety.vtk");
		series3DWriter->Update();


		using LabelGeometryImageFilterType = itk::LabelGeometryImageFilter<Image3DType>;
		LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
		labelGeometryImageFilter->SetInput(connComp1->GetOutput());
		//labelGeometryImageFilter->SetIntensityInput(connComp->GetOutput());

		// These generate optional outputs.
		labelGeometryImageFilter->CalculatePixelIndicesOn();
		labelGeometryImageFilter->CalculateOrientedBoundingBoxOn();

		labelGeometryImageFilter->Update();
		LabelGeometryImageFilterType::LabelsType labelsAll= labelGeometryImageFilter->GetLabels();
		
		std::cout << labelsAll.size() << std::endl;
		std::cout << labelGeometryImageFilter->GetCentroid(labelsAll.size()-1) << std::endl;
		
		LabelGeometryImageFilterType::LabelPointType startPoint = labelGeometryImageFilter->GetCentroid(labelsAll.size() - 1);
		//LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
		
		std::cout << startPoint[0] << std::endl;
		std::cout << startPoint[1] << std::endl;
		std::cout << startPoint[2] << std::endl;
	//===============================================================================================================

		using ConnectedFilterType =itk::ConfidenceConnectedImageFilter<Image3DType, Image3DType>;
		ConnectedFilterType::Pointer confidenceConnected = ConnectedFilterType::New();
		confidenceConnected->SetInput(image3D); // daæ tu obraz po mno¿eniu zeby bylo przyciête
		//confidenceConnected->SetMultiplier(0.57);  //Multiplier = 0.57, iteracje 1
		//confidenceConnected->SetNumberOfIterations(1);
		confidenceConnected->SetMultiplier(1.3);  //Multiplier = 0.57, iteracje 1
		confidenceConnected->SetNumberOfIterations(1);
		confidenceConnected->SetReplaceValue(99);
		confidenceConnected->SetInitialNeighborhoodRadius(2);
		ConnectedFilterType::IndexType index;
		index[0] = (int)(startPoint[0]) + 2; // 152;//129;//152;
		index[1] = (int)(startPoint[1])+2; //168;//129;// 168;
		index[2] = (int)(startPoint[2])+2;//48;// 48;

		std::cout << "Indeksy:" << std::endl;
		std::cout << image3D->GetPixel(index) << std::endl;
		std::cout << index[0] << std::endl;
		std::cout << index[1] << std::endl;
		std::cout << index[2] << std::endl;
		std::cout << "KONIEC" << std::endl;

		//index[0] = 111;//117;//129;//152;
		//index[1] = 130;// 158;//129;// 168;
		//index[2] = 47;// 48;
		//ConnectedFilterType::IndexType index2;
		//index2[0] = 111;//129;//152;
		//index2[1] = 168;//129;// 168;
		//index2[2] = 48;// 48;
	confidenceConnected->AddSeed(index);
	//confidenceConnected->AddSeed(index2);
	confidenceConnected->Update();
	series3DWriter->SetInput(confidenceConnected->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9_img3D_rozrost-thresh.vtk");
	series3DWriter->Update();
	auto pixelVal=image3D->GetPixel(index);
	
	std::cout << image3D_mnozenie->GetPixel(index) << std::endl;

	//multiply->SetInput1(confidenceConnected->GetOutput());
	//multiply->SetInput2(image3D_mnozenie);
	////image3D_mnozenie = multiply->GetOutput();;
	//multiply->Update();
	//series3DWriter->SetInput(multiply->GetOutput());
	//series3DWriter->SetFileName("..\\wyniki\\img3D_mnozenie_v2.vtk");
	//series3DWriter->Update();


	using InvertIntensityImageFilterType = itk::InvertIntensityImageFilter<Image3DType>;

	InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
	invertIntensityFilter->SetInput(confidenceConnected->GetOutput());
	invertIntensityFilter->SetMaximum(1);
	invertIntensityFilter->Update();
	series3DWriter->SetInput(invertIntensityFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9a_img3D_invert.vtk");
	series3DWriter->Update();

	using WindowingImageFilter = itk::IntensityWindowingImageFilter<Image3DType>;
	WindowingImageFilter::Pointer windowingImageFilter = WindowingImageFilter::New();
	windowingImageFilter->SetInput(invertIntensityFilter->GetOutput());
	windowingImageFilter->SetWindowMaximum(1);
	windowingImageFilter->SetWindowMinimum(-1);
	windowingImageFilter->SetOutputMinimum(0);
	windowingImageFilter->SetOutputMaximum(1);
	windowingImageFilter->Update();
	series3DWriter->SetInput(windowingImageFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9b_img3D_window.vtk");
	series3DWriter->Update();

	

	Image3DType::Pointer windowImage = Image3DType::New();
	Image3DType::Pointer erode = Image3DType::New();

	

	reader3D->SetFileName("..\\wyniki\\9b_img3D_window.vtk");
	reader3D->Update();
	windowImage = reader3D->GetOutput();

	reader3Db->SetFileName("..\\wyniki\\4_img3D_erozja.vtk");
	reader3Db->Update();
	erode = reader3Db->GetOutput();

	index[0] = 111;//117;//129;//152;
	index[1] = 130;// 158;//129;// 168;
	index[2] = 47;// 48;

	std::cout << windowingImageFilter->GetOutput()->GetPixel(index) << std::endl;
	std::cout << reader3Db->GetOutput()->GetPixel(index) << std::endl;

	multiply->SetInput1(windowImage);
	multiply->SetInput2(erode);
	//image3D_mnozenie = multiply->GetOutput();;
	multiply->Update();
	series3DWriter->SetInput(multiply->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9c_img3D_mnozenie_komory.vtk");
	series3DWriter->Update();



	using ConnCompFilterType = itk::ConnectedComponentImageFilter<Image3DType, Image3DType>;

	ConnCompFilterType::Pointer connComp = ConnCompFilterType::New();
	connComp->SetInput(multiply->GetOutput());
	connComp->SetBackgroundValue(0);
	connComp->Update();

	series3DWriter->SetInput(connComp->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9d_img3D_label.vtk");
	series3DWriter->Update();


	using RelabelComponentFilterType = itk::RelabelComponentImageFilter<Image3DType, Image3DType>;
	RelabelComponentFilterType::Pointer relabel = RelabelComponentFilterType::New();
	relabel->SetInput(connComp->GetOutput());
	relabel->SetMinimumObjectSize(100);
	
	relabel->Update();

	series3DWriter->SetInput(relabel->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9e_img3D_relabel.vtk");
	series3DWriter->Update();


	using FilterType = itk::BinaryThresholdImageFilter<Image3DType, Image3DType>;
	FilterType::Pointer thresholder2 = FilterType::New();
	thresholder2->SetInput(relabel->GetOutput());
	thresholder2->SetLowerThreshold(2);
	thresholder2->SetUpperThreshold(2);
	thresholder2->Update();

	series3DWriter->SetInput(thresholder2->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9f_img3D_binaryzacja1komora.vtk");
	series3DWriter->Update();
	
	//flip
	using FlipImageFilterType = itk::FlipImageFilter<Image3DType>;

	FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
	flipFilter->SetInput(thresholder2->GetOutput());


	FlipImageFilterType::FlipAxesArrayType flipAxes;

	flipAxes[0] = true;
	flipAxes[1] = false;
	flipAxes[2] = false;
	flipFilter->SetFlipAxes(flipAxes);
	
	flipFilter->Update();

	series3DWriter->SetInput(flipFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9g_img3D_flip.vtk");
	series3DWriter->Update();

	//add
	//reader3Dc->SetFileName("..\\wyniki\\img3D_binaryzacja1komora.vtk");
	////reader3Dc->SetFileName("..\\wyniki\\img3D.vtk");
	//reader3Dc->Update();
	//Image3DType::Pointer komora13D = Image3DType::New();
	//komora13D = reader3Dc->GetOutput();

	//reader3Dd->SetFileName("..\\wyniki\\img3D_flip.vtk");
	//reader3Dd->Update();
	Image3DType::Pointer flipImage = Image3DType::New();
	//flip = reader3Dd->GetOutput();

	flipImage = flipFilter->GetOutput();
	flipImage->CopyInformation(thresholder2->GetOutput());
	

	
	




	using ResampleFilterType = itk::ResampleImageFilter<Image3DType, Image3DType>;
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();

	using TransformType = itk::AffineTransform<double, 3>;
	TransformType::Pointer transform = TransformType::New();
	TransformType::OutputVectorType translation;
	translation[0] = -8; // X translation in millimeters
	translation[1] = 0; // Y translation in millimeters
	translation[2] = 0; // Z translation in millimeters
	transform->Translate(translation);

	resampler->SetInput(flipImage);
	resampler->SetTransform(transform);
	
	resampler->SetSize(flipImage->GetLargestPossibleRegion().GetSize());
	resampler->SetOutputOrigin(flipImage->GetOrigin());
	resampler->SetOutputSpacing(flipImage->GetSpacing());
	resampler->SetOutputDirection(flipImage->GetDirection());
	resampler->SetDefaultPixelValue(100);

	series3DWriter->SetInput(resampler->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9j_img3D_resample.vtk");
	series3DWriter->Update();




	using AddImageFilterType = itk::AddImageFilter<Image3DType, Image3DType>;

	AddImageFilterType::Pointer addFilter = AddImageFilterType::New();


	Image3DType::Pointer flippedImage = Image3DType::New();
	flippedImage = resampler->GetOutput();
	flippedImage->CopyInformation(thresholder2->GetOutput());

	addFilter->SetInput1(flippedImage);
	addFilter->SetInput2(thresholder2->GetOutput());


	series3DWriter->SetInput(addFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9h_img3D_dodawanie.vtk");
	series3DWriter->Update();



	WindowingImageFilter::Pointer windowingFilter = WindowingImageFilter::New();
	windowingFilter->SetInput(addFilter->GetOutput());
	windowingFilter->SetWindowMaximum(1000);
	windowingFilter->SetWindowMinimum(0);
	windowingFilter->SetOutputMinimum(0);
	windowingFilter->SetOutputMaximum(1000);
	windowingFilter->Update();
	series3DWriter->SetInput(windowingFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9i_img3D_dodawanie_window.vtk");
	series3DWriter->Update();



	Image3DType::Pointer maskImage = Image3DType::New();
	maskImage = windowingFilter->GetOutput();
	maskImage->CopyInformation(image3D);

	AddImageFilterType::Pointer addFilter2 = AddImageFilterType::New();
	addFilter2->SetInput1(maskImage);
	addFilter2->SetInput2(image3D);


	series3DWriter->SetInput(addFilter2->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9j_img3D_nalozone_maski.vtk");
	series3DWriter->Update();
	

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
