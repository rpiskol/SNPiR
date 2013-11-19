import java.util.ArrayList;
import java.util.Arrays;




public class convertCoordinates {

	private static String convertLine(String line) {
		// TODO Auto-generated method stub
		//System.out.println(line);
		String[] fields = line.split("\t");
		//System.out.println(fields[2]);
		String[] chrom = fields[2].split("[-]");
		//System.out.println(Arrays.toString(chrom));
		//System.out.println(line.charAt(1));
		if(line.charAt(0) == '@'){
			return line;
		}
		else if ( fields[5].equals("*") && ( fields[2].length() != 1 || Integer.parseInt(fields[3]) != 0 ) ){
			//	System.out.println("malformed read");
				//return (fields[0] + "malformed read");
				//System.exit(0);
				fields[2] = "*";
				fields[3] = "0";
				String newline = Arrays.toString(fields).replace(", ", "\t").replaceAll("[\\[\\]]", "");
				return newline;
				
			}
		//else if (chrom.length == 7 && )
		else if (chrom.length == 1) {
			return line;
		}
		//if(chrom.length != 6){
		//	return line;
		//}
		else if(chrom.length == 100){
			
			int start1 = Integer.parseInt(chrom[1]);
			int stop1 = Integer.parseInt(chrom[2]);
			int start2 = Integer.parseInt(chrom[3]);
			int beginning = Integer.parseInt(fields[3]);
			int read1start = start1 + beginning - 1;
			//System.out.println(read1start);
            int read1length = stop1 - read1start + 1;
            int Nlength = start2-stop1-1;
            String cigar = fields[5];
            String newcigar = "";
            
            
            //System.out.println(cigarletters[1]);
            
            if(cigar.equals("*")){
            	newcigar = "*";
            }
            else if (read1length <= 0) {
				newcigar = cigar;
				//read1start = start2 - (read1length + 1);
				read1start = start2 + beginning - stop1 + start1 - 2;
				//System.out.println(read1start);
			}
            else{
            	String[] cigarnumsTmp = cigar.split("[MIDNSHP]");
                int[] cigarnums = new int[cigarnumsTmp.length];
                for(int i = 0; i < cigarnumsTmp.length; i++){
                	cigarnums[i] = Integer.parseInt(cigarnumsTmp[i]);
                }
                String cigarlettersTmp[] = cigar.split("[0-9]+");//.toString().toCharArray();
                String cigarletters[] = Arrays.copyOfRange(cigarlettersTmp, 1, cigarlettersTmp.length);
            	int currentpos = 0;
            	boolean found = false;
            	
            	for(int i = 0; i < cigarnums.length; i++){
            		if(!cigarletters[i].equals("I") && !cigarletters[i].equals("S") && !cigarletters[i].equals("H")){
            			//System.out.println("adding cignum" + cigarnums[i] + " for " + );
            			currentpos += cigarnums[i];
            		}
            		
            		if(currentpos > read1length && !found){
            			found = true;
            			int len1 = cigarnums[i] - (currentpos - read1length);
            			int len2 = cigarnums[i] - len1;
            			newcigar = newcigar + String.valueOf(len1) + cigarletters[i] + String.valueOf(Nlength) + "N" + String.valueOf(len2 + cigarletters[i]);
            		}
            		else{
            			newcigar = newcigar + String.valueOf(cigarnums[i]) + cigarletters[i];
            		}
            	}
            }                
            String sep = "\t";
            String newLine = fields[0] + sep + fields[1] + sep + chrom[0] + sep + read1start + sep + fields[4] + sep + newcigar + sep;
            for(int i = 6; i < fields.length-1; i++){
            	newLine += fields[i]+sep;
            }
            newLine += fields[fields.length-1];
            
            return newLine;

		}
		else{
			
			String newline = parseMultiexon(fields,chrom);
			return newline;
			
			
		}
	}

	
	
	private static String parseMultiexon(String[] fields, String[] chrom) {
		ArrayList<Integer> exonStarts = new ArrayList<Integer>();
		ArrayList<Integer> exonEnds = new ArrayList<Integer>();
		ArrayList<Integer> cmLength = new ArrayList<Integer>();
		//ArrayList<Integer> exonLength = new ArrayList<Integer>();
		
		int readStartRel = Integer.parseInt(fields[3]);
		String newcigar = "";
		
		/*
		 * get numbers and characters from cigar string
		 */
		String cigar = fields[5];
		String[] cigarnumsTmp = cigar.split("[MIDNSHP]");
        int[] cigarnums = new int[cigarnumsTmp.length];
        int cigarsum = 0;
        for(int i = 0; i < cigarnumsTmp.length; i++){
        	cigarnums[i] = Integer.parseInt(cigarnumsTmp[i]);
        	cigarsum += cigarnums[i]; 
        }
        String cigarlettersTmp[] = cigar.split("[0-9]+");//.toString().toCharArray();
        String cigarletters[] = Arrays.copyOfRange(cigarlettersTmp, 1, cigarlettersTmp.length);
		
        /*
		 * get exon starts and ends as well as individual and cummulative exon lengths
		 */
		int len = 0;
		for(int i = 1; i < chrom.length; i++){
			//System.err.println(i + "\t" + Arrays.toString(chrom));
			if( i % 2 == 0){
				exonEnds.add(Integer.parseInt(chrom[i]));
				len += exonEnds.get(exonEnds.size()-1) - exonStarts.get(exonStarts.size()-1) +1;
				//exonLength.add(new Integer(exonEnds.get(exonEnds.size()-1) - exonStarts.get(exonStarts.size()-1) +1));
				cmLength.add(new Integer(len));
			}
			else{
				exonStarts.add(Integer.parseInt(chrom[i]));
			}
		}
		//System.err.println(cmLength);
		/*
		 * find the exon that contains the read start
		 */
		int readStartIdx = -1;
		int readStartAbs = -1;
		if(readStartRel <= cmLength.get(0)){
			readStartIdx = 0;
			readStartAbs = exonStarts.get(0) + readStartRel - 1;
		}
		for(int i = 1; i < cmLength.size(); i++){
			//System.err.println(readStartRel + "\t" +  cmLength.get(i)  + "\t" +  readStartRel  + "\t" +   (cmLength.get(i)-cmLength.get(i-1)));
			if(readStartRel <= cmLength.get(i) && readStartRel > cmLength.get(i-1)){
				readStartIdx = i;
				readStartAbs = exonStarts.get(i) + (readStartRel-cmLength.get(i-1))-1;
			}
		}

        
        //System.err.println(readStartIdx);
        if(cigar.equals("*")){
        	newcigar = "*";
        }
        else{
        	int currCigarPos = 0;
        	int currReadPos = exonEnds.get(readStartIdx) - readStartAbs + 1;
            int readRemain = cigarsum;
        	boolean found = false;

        	for(int i = 0; i < cigarnums.length; i++){
        		
        		if(!cigarletters[i].equals("I") && !cigarletters[i].equals("S") && !cigarletters[i].equals("H")){
        			//System.out.println("adding cignum" + cigarnums[i] + " for " + );
        			currCigarPos += cigarnums[i];
        			//readRemain -= cigarnums[i];
        		}
        		
        		while(currCigarPos > currReadPos && !found && readStartIdx < exonEnds.size()-1){

        			//found = true;
        			int len1 = cigarnums[i] - (currCigarPos - currReadPos);
        			cigarnums[i] -= len1;
        			int Nlength = exonStarts.get(readStartIdx+1) - exonEnds.get(readStartIdx) - 1 ;
        			
        			
        			if(len1>0){//only include current cigarchar and length if it is larger than zero
        				newcigar = newcigar + String.valueOf(len1) + cigarletters[i] + String.valueOf(Nlength) + "N";// + String.valueOf(len2 + cigarletters[i]);
        			}
        			else{
        				newcigar = newcigar + String.valueOf(Nlength) + "N";// + String.valueOf(len2 + cigarletters[i]);
        			}
        			readRemain -= len1;
        			//System.err.println(readStartIdx + "\t" + readRemain);
        			readStartIdx++;
        			//System.err.println(exonEnds.get(readStartIdx) - exonStarts.get(readStartIdx) + 1);
        			if(exonEnds.get(readStartIdx) - exonStarts.get(readStartIdx) + 1 < readRemain){
        				currReadPos += (exonEnds.get(readStartIdx) - exonStarts.get(readStartIdx) + 1);
        				//System.err.println("here");
        			}
        			else{
        				currReadPos += readRemain;
        			}
        			
        		}
        		//else{
        		if(cigarnums[i]>0){
        			newcigar = newcigar + String.valueOf(cigarnums[i]) + cigarletters[i];
        		}
        		//}
        	}
        }
//		System.err.println(readStartRel);
//		System.err.println(readStartIdx);
//		System.err.println(readStartAbs);
		//System.err.println(exonStarts);
		//System.err.println(exonEnds);
//		System.err.println(exonLength);
//		System.err.println(cmLength);
//		System.err.println(cigar);
//		System.err.println(newcigar);
		
		String sep = "\t";
		String newLine = fields[0] + sep + fields[1] + sep + chrom[0] + sep + readStartAbs + sep + fields[4] + sep + newcigar + sep;
        for(int i = 6; i < fields.length-1; i++){
        	newLine += fields[i]+sep;
        }
        newLine += fields[fields.length-1];
        
		return newLine;
	}



	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
	       java.io.BufferedReader stdin = new java.io.BufferedReader(new java.io.InputStreamReader(System.in));
	       String line;
	       
	       
			while((line = stdin.readLine()) !=null){
				//System.out.println();
				 //System.err.println('>'+line);
				 String convertedLine = convertLine(line);
				 
				 System.out.println(convertedLine);
			}
		}	  
		catch (java.io.IOException e) { System.out.println(e); }
		catch (NumberFormatException e) { System.out.println(e); }


	}


}
