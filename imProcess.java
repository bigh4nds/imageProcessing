import com.jhlabs.image.*;
import java.awt.*;
import java.awt.image.*;
import javax.imageio.*;
import java.io.*;
import java.util.*;

class imProcess {

	public static void printUsage() {
		System.out.println("Usage: java imProcess -input <filename> [options]");
		System.out.println("-output <filename>");
		System.out.println("-brightness <float>");
		System.out.println("-edgedetect");
		System.out.println("-blur <float>");

		System.out.println("-contrast <float>");
		System.out.println("-saturation <float>");
		System.out.println("-sharpen <float>");
		System.out.println("-randomdither");
		System.out.println("-ordereddither");
		System.out.println("-mosaic <imagefolder>");
		System.exit(1);
	}

	public static void main(String[] args) {
		int i = 0;
		BufferedImage src = null, dst = null;
		BufferedImage tmp = null;	// used for swapping src and dst buffer
		int width, height;			// image width, height

		String arg;
		String outputfilename = "output.png";		// default output filename
		
		if (args.length < 2) {
			printUsage();
		}
		
		// parse command line options, and call approrpiate member functions
		while (i < args.length && args[i].startsWith("-")) {
			arg = args[i++];

			if (arg.equals("-input")) {

				String inputfile = args[i++];
				try {
					src = ImageIO.read(new File(inputfile));
				} catch (IOException e) {
				}
				width = src.getWidth();
				height = src.getHeight();
				dst = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
				continue;

			} else if (arg.equals("-output")) {

				outputfilename = args[i++];
				System.out.println("Output file: " + outputfilename);
				continue;

			} else if (arg.equals("-brightness")) {

				float brightness = Float.parseFloat(args[i++]);
				System.out.println("Set brightness: " + brightness);
				Brighten(src, dst, brightness);

			} else if (arg.equals("-contrast")) {

				float contrast = Float.parseFloat(args[i++]);
				System.out.println("Set contrast: " + contrast);

				AdjustContrast(src, dst, contrast);

			} else if (arg.equals("-saturation")) {

				float saturation = Float.parseFloat(args[i++]);
				System.out.println("Set saturation: " + saturation);

				AdjustSaturation(src, dst, saturation);

			} else if (arg.equals("-randomdither")) {

				System.out.println("Generated random dithering");

				RandomDither(src, dst);

			} else if (arg.equals("-ordereddither")) {

				System.out.println("Generate ordered dithering");

				OrderedDither(src, dst);
			} else if (arg.equals("-blur")) {

				float radius = Float.parseFloat(args[i++]);
				System.out.println("Set blur radius: " + radius);
				
				Blur(src, dst, radius);

			} else if (arg.equals("-sharpen")) {

				float sharpness = Float.parseFloat(args[i++]);
				System.out.println("Set sharpness : " + sharpness);

				Sharpen(src, dst, sharpness);

			} else if (arg.equals("-edgedetect")) {

				System.out.println("Apply edge detector");

				EdgeDetect(src, dst);

			} else if (arg.equals("-mosaic")) {

				String mosaicfolder = args[i++];
				System.out.println("Base image folder: " + mosaicfolder);

				Mosaic(src, dst, mosaicfolder);

			} else {
				printUsage();
			}
			// swap src and dst to prepare for the next operation
			tmp = src; src = dst; dst = tmp;
		}
		if (i != args.length) {
			System.out.println("there are unused arguments");
		}
		// write the output image to disk file
		File outfile = new File(outputfilename);
		try {
			ImageIO.write(src, "png", outfile);
		} catch(IOException e) {
		}
	}

	// Change the brightness of an image
	// brightness is a scaling factor 
	// Use this function as an example. There is nothing you need to change here
	public static void Brighten(BufferedImage src, BufferedImage dst, float brightness) {

		int width = src.getWidth();
		int height = src.getHeight();

		// a buffer that stores the destination image pixels
		int[] pixels = new int[width * height];
	
		// get the pixels of the source image	
		src.getRGB(0, 0, width, height, pixels, 0, width);

		int i;
		int a, r, g, b;
		for(i = 0; i < width * height; i ++) {
			Color rgb = new Color(pixels[i]);
			// a color is represented as an integer (4 bytes); 
			// each of the alpha, red, green, blue channels is stored in one byte in order;
			// you can use the Color class to extract the value of each individual channel
			// or composite a new integer color from the separated components
			a = rgb.getAlpha();
			r = rgb.getRed();
			g = rgb.getGreen();
			b = rgb.getBlue();
			r = PixelUtils.clamp((int)((float)r * brightness));
			g = PixelUtils.clamp((int)((float)g * brightness));
			b = PixelUtils.clamp((int)((float)b * brightness));

			pixels[i] = new Color(r, g, b, a).getRGB();
		}

		// write pixel values to the destination image
		dst.setRGB(0, 0, width, height, pixels, 0, width);

	}

	// change the contrast of an image
	// contrast = 0 gives a medium gray (0.5, 0.5, 0.5) image
	// constrat = 1 gives the original image
	public static void AdjustContrast(BufferedImage src, BufferedImage dst, float contrast) {

		int width = src.getWidth();
		int height = src.getHeight();

		// a buffer that stores the destination image pixels
		int[] pixels = new int[width * height];
	
		// get the pixels of the source image	
		src.getRGB(0, 0, width, height, pixels, 0, width);

		int i;
		int avgGray = 256/2;
		int a, r, g, b;
		for(i = 0; i < width * height; i ++) {
			Color rgb = new Color(pixels[i]);
			a = rgb.getAlpha();
			r = rgb.getRed();
			g = rgb.getGreen();
			b = rgb.getBlue();
			r = PixelUtils.interpolate(avgGray, r, contrast);
			g = PixelUtils.interpolate(avgGray, g, contrast);
			b = PixelUtils.interpolate(avgGray, b, contrast);
			pixels[i] = new Color(r, g, b, a).getRGB();
		}

		// write pixel values to the destination image
		dst.setRGB(0, 0, width, height, pixels, 0, width);
	}

	// change the saturation of an image
	// saturation = 0 gives a gray scale version of the image
	// saturation = 1 gives the original image
	public static void AdjustSaturation(BufferedImage src, BufferedImage dst, float saturation) {

		int width = src.getWidth();
		int height = src.getHeight();

		// a buffer that stores the destination image pixels
		int[] pixels = new int[width * height];
	
		// get the pixels of the source image	
		src.getRGB(0, 0, width, height, pixels, 0, width);

		int i;
		double L;
		int a, r, g, b;
		for(i = 0; i < width * height; i ++) {
			Color rgb = new Color(pixels[i]);
			a = rgb.getAlpha();
			r = rgb.getRed();
			g = rgb.getGreen();
			b = rgb.getBlue();
			L = 0.30*r + 0.59*g + 0.11*b;
			r = PixelUtils.interpolate((int)L, r, saturation);
			g = PixelUtils.interpolate((int)L, g, saturation);
			b = PixelUtils.interpolate((int)L, b, saturation);
			pixels[i] = new Color(r, g, b, a).getRGB();
		}

		// write pixel values to the destination image
		dst.setRGB(0, 0, width, height, pixels, 0, width);
	}

	// blur an image
	// use the GaussianFilter from jhlabs to perform Gaussian blur
	// This function is a given, and there is nothing you need to change here
	public static void Blur(BufferedImage src, BufferedImage dst, float radius) {

		GaussianFilter filter = new GaussianFilter();
		filter.setRadius(radius);
		filter.filter(src, dst);
	}

	// sharpen an image
	// sharpness sets the amount of sharpening
	// use the ConvolveFilter from jhlabs and the sharpening matrix we covered in class to perform this operation
	public static void Sharpen(BufferedImage src, BufferedImage dst, float sharpness) {

		ConvolveFilter filter = new ConvolveFilter();
		float[] matrix =      {0f, -sharpness, 0f,
		  		        -sharpness, 1+4*sharpness, -sharpness,
		  		              0f,  -sharpness,  0f};
		Kernel kernel = new Kernel(3, 3, matrix);
		filter.setKernel(kernel);
		filter.filter(src, dst);	
	}

	// detect edge features of an image
	// use the EdgeFilter from jhlabs
	// This function is a given, and there is nothing you need to change here
	public static void EdgeDetect(BufferedImage src, BufferedImage dst) {

		EdgeFilter filter = new EdgeFilter();
		filter.filter(src, dst);
	}

	// random dithering
	// compare each image pixel against a random threshold to quantize it to 0 or 1
	// ignore the color, and just use the luminance of a pixel to do the dithering
	// your output should be a binary (black-white) image
	public static void RandomDither(BufferedImage src, BufferedImage dst) {
		
		int width = src.getWidth();
		int height = src.getHeight();

		// a buffer that stores the destination image pixels
		int[] pixels = new int[width * height];
	
		// get the pixels of the source image	
		src.getRGB(0, 0, width, height, pixels, 0, width);

		int i;
		double L;
		int a, r, g, b;
		double threshold;
		for(i = 0; i < width * height; i ++) {
			threshold = Math.random();
			Color rgb = new Color(pixels[i]);
			a = rgb.getAlpha();
			r = rgb.getRed();
			g = rgb.getGreen();
			b = rgb.getBlue();
			L = (0.30*r + 0.59*g + 0.11*b)/256;
			if(L > threshold)
				pixels[i] = new Color(255, 255, 255, a).getRGB();
			else
				pixels[i] = new Color(0, 0, 0, a).getRGB();		
		}

		// write pixel values to the destination image
		dst.setRGB(0, 0, width, height, pixels, 0, width);
	}
	
	// ordered dithering
	// compare each image pixel against a pseudo-random threshold to quantize it to 0 or 1
	// in this case, the pseudo random number is given by a 4x4 Bayers matrix
	// ignore the color, and just use the luminance of a pixel to do the dithering
	// your output should be a binary (black-white) image
	public static void OrderedDither(BufferedImage src, BufferedImage dst) {
		
		final float[][] Bayers = {{15/16.f,  7/16.f,  13/16.f,   5/16.f},
								  {3/16.f,  11/16.f,   1/16.f,   9/16.f},
								  {12/16.f,  4/16.f,  14/16.f,   6/16.f},
								  { 0,      8/16.f,    2/16.f,  10/16.f}};

		int width = src.getWidth();
		int height = src.getHeight();

		// a buffer that stores the destination image pixels
		int[] pixels = new int[width * height];
	
		// get the pixels of the source image	
		src.getRGB(0, 0, width, height, pixels, 0, width);

		int i;
		int y = 0;
		double L;
		int a, r, g, b;
		double threshold;
		for(i = 0; i < width * height; i ++) {
			threshold = Bayers[y%4][i%4];
			if((i%(width-1) == 0))
				y++;
			Color rgb = new Color(pixels[i]);
			a = rgb.getAlpha();
			r = rgb.getRed();
			g = rgb.getGreen();
			b = rgb.getBlue();
			L = (0.30*r + 0.59*g + 0.11*b)/256;
			if(L > threshold)
				pixels[i] = new Color(255, 255, 255, a).getRGB();
			else
				pixels[i] = new Color(0, 0, 0, a).getRGB();	
		}

		// write pixel values to the destination image
		dst.setRGB(0, 0, width, height, pixels, 0, width);
	}

	// generate image Mosaics
	// mosaicfolder specifies a subfolder containing a collection of images
	// to be used for producing Mosaic
	public static void Mosaic(BufferedImage src, BufferedImage dst, String mosaicfolder) {

		int width = src.getWidth();
		int height = src.getHeight();

		// load all mosaic images from the specified subfolder
		File folder = new File(mosaicfolder);
		File files[] = folder.listFiles();

		int i;
		int w = 0, h = 0;
		int num = files.length;

		// mpixels stores the pixels of each mosaic image read from a disk file
		int[][] mpixels = new int[num][];

		for (i = 0; i < files.length; i++) {
			if (!files[i].isFile()) continue;
			BufferedImage mosaic = null;
			try {
				mosaic = ImageIO.read(files[i]);
			} catch (IOException e) {
			}
			if (w == 0) {
				w = mosaic.getWidth();
				h = mosaic.getHeight();
			} else {
				if (mosaic.getWidth() != w || mosaic.getHeight() != h) {
					System.out.println("mosaic images must be of the same size.");
					System.exit(1);
				}
			}
			mpixels[i] = new int[w*h];

			// get pixels from the buffered image
			mosaic.getRGB(0, 0, w, h, mpixels[i], 0, w);
		}
		System.out.println("" + num + " mosaic images (" + w + "," + h + ") loaded.");

		/* === YOUR WORK HERE === */
		int[] srcPixels = new int[width*height];
		float[] Nr = new float[num]; // Precalculated
		float[] Ng = new float[num]; // denominators
		float[] Nb = new float[num]; // for possible matches
		int bestMatch = 0; // bestMatch mpixels index
		int W = w; // Used to handle
		int H = h; // cutoff at edges
		float tffsq = 255*255; // Used as divisor to scale d's to value between 0.0 & 1.0	
		for(i=0;i<num;i++) // Loop through all mosaic images
		{
			for(int k=0;k<w*h;k++) // Loop though all pixels of current mosaic image
			{
				Color mos = new Color(mpixels[i][k]); 
				Nr[i] += (mos.getRed()*mos.getRed())/tffsq;     // Precomputing
				Ng[i] += (mos.getGreen()*mos.getGreen())/tffsq; // denominator (Sum(Mk^2))
				Nb[i] += (mos.getBlue()*mos.getBlue())/tffsq;   // for A & d equations
			}
		}
		for(int y=0;y<height;y=y+h) // Go through height of srcImage
		{
			if(y+h > height)        // Check for uneven patterning and cut off to fit
				H = height - (y+h);
			W = w;
			for(int x=0;x<width;x=x+w) // Go through width of srcImage
			{
				if(x+w > width)        // Check for uneven patterning and cut off to fit
					W = width - (x+w);
				float lowest_d = Float.MAX_VALUE;      // Keep track of lowest d value
				float Ar=0, Ag=0, Ab=0; // bestMatch alphas			
				for(i=0;i<num;i++) // Check through all possible images
				{
					src.getRGB(x, y, W, H, srcPixels, 0, W); // Get srcImage's pixels at current WxH square
					float dr=0, dg=0, db=0, d=0; // Reinitialize d's to 0
					for(int k=0;k<W*H;k++) { // Loop through all pixels in WxH square
						Color src_rgb = new Color(srcPixels[k]);  
						Color mos_rgb = new Color(mpixels[i][k]); 
						int mos_r = mos_rgb.getRed();
						int mos_g = mos_rgb.getGreen();
						int mos_b = mos_rgb.getBlue();
						int src_r = src_rgb.getRed();
						int src_g = src_rgb.getGreen();
						int src_b = src_rgb.getBlue();
						dr += (src_r*mos_r)/tffsq; // d
						dg += (src_g*mos_g)/tffsq; // Numerator
						db += (src_b*mos_b)/tffsq; // Summations
					}
					d = (((-(dr*dr)/Nr[i]) + (-(dg*dg)/Ng[i]) + (-(db*db)/Nb[i])) * (float)(Math.random()+1)); // negation, square, and division by Sum([Mk[i] pixels)^2
					if(d < lowest_d)  // Update current bestMatch
					{
						lowest_d = d;
						Ar = dr/Nr[i];
						Ag = dg/Ng[i];
						Ab = db/Nb[i];
						bestMatch = i;
					}
				}
				for(int j=0;j<w*h;j++)
				{
					Color mos_rgb = new Color(mpixels[bestMatch][j]);
					float mos_r = mos_rgb.getRed();
					float mos_g = mos_rgb.getGreen();
					float mos_b = mos_rgb.getBlue();
					mos_r = (float)PixelUtils.clamp((int)(Ar*(mos_r)))/255; // scaling
					mos_g = (float)PixelUtils.clamp((int)(Ag*(mos_g)))/255; // RGB values
					mos_b = (float)PixelUtils.clamp((int)(Ab*(mos_b)))/255; // by alpha	
					srcPixels[j] = new Color(mos_r, mos_g, mos_b, mos_rgb.getAlpha()/255).getRGB(); // update srcPixel with chosen & scaled RGB values
				}
				dst.setRGB(x, y, W, H, srcPixels, 0, W); // Update dstImage with new WxH pixel block
			}
		}
	}
}
