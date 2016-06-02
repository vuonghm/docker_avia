#!/usr/bin/perl
use strict;
use lib '/users/abcc/vuonghm/scripts.dir/perl_modules/JSON-2.90/lib/';
use Data::Dumper;
use JSON;
#read in .STATS file to build original table

#read in uniprot.json file to build secondary table
if (-e "uniprot.json"){
	open (FILE,"<uniprot.json") or die "Cannot open uniprot.json\n";
	my $uni_txt = <FILE>;
	my $uniprot_decoded=json_decode($uni_txt);
	print $uniprot_decoded;
# 	$arr=json_decode($uni_json,true);
# 	$uniprot_header=array();
# 	$table_body='';
# 	$patterns=array('/\(/','/\)/');
# 	foreach ($arr as $mutation =>$obj){
# 		foreach ($obj as $key => $element){
# 			if (!array_key_exists("$key", $uniprot_header)) {
# 				$uniprot_header[$key]=1;
# 			}
# 		}
# 	}
# 	ksort($arr);
# 	foreach ($arr as $mutation =>$obj){
# 		if (preg_match("/^\w+:([A-Z]*)(\d+)([A-Z]*)/",$mutation,$genes)){
# 			$gene=preg_replace("/:/","",$genes[0]);
# 			$table_body.="<tr><td valign=\"top\">$mutation</td>";
# 			// if (in_array("$mutation", $jmol_arr)){
# 			// 	$table_body.="<td><a href=\"/apps/site/jmol/?id=$id\"><img src=\"/modules/resource_manager/library/images/icon_link_out.png\"/></a></td>";
# 			// }else{
# 			// 	$table_body.="<td>&nbsp;</td>";
# 			// }
# 		}else{
# 			$table_body.="<tr><td valign=\"top\">$mutation</td>";//<td>&nbsp;</td>";
# 		}
# 		foreach ($uniprot_header as $key => $element){
# 			if (!array_key_exists("$key", $obj)) {
# 				$table_body.= "<td>&nbsp;</td>
# 				";
# 			}else{
# 				if (preg_match("/accession/",$key)){
# 					// $accessions=implode(",",$arr[$mutation][$key]);
# 					$table_body.="<td style=\"min-width:20px;\" valign=\"top\">";$addon1='
# 					<td valign="top">';
# 					foreach ($arr[$mutation][$key] as $keypair=>$acc){
# 						$table_body.="<a title=\"Link out to Uniprot $acc\" href=\"http://www.uniprot.org/uniprot/$acc\" target=\"_blank\" style=\"text-decoration:underline;color: #006DB5;\">$acc</a>, ";
# 						$addon1.="
# 						<a href=\"http://www.proteinmodelportal.org/query/uniprot/$acc\" title=\"Link to protein model portal \" target=\"_blank\" style=\"text-decoration:underline;color: #006DB5;\">$acc</a>, 
# 						";
# 					}
# 					$addon1.="</td>";
# 					$addon2="<td valign=\"top\" align=\"center\"><a href=\"http://reactome.org/cgi-bin/link?SOURCE=Uniprot&amp;ID=$acc\" title=\"Link out to Reactome for $acc\" target=\"_blank\" style=\"text-decoration:underline;color: #006DB5;\"><img src=\"/modules/resource_manager/library/images/icon_link_out.png\"/></a></td>";
# 					$addon3='<td>&nbsp;</td>';
# 					$table_body.="</td>$addon1$addon2
# 					";
# 				}else{
# 					$table_body.="<td valign=\"top\">".preg_replace($patterns,'',implode(',',$arr[$mutation][$key]))."</td>
# 					";
# 				}
# 			}
# 		}
# 		$table_body.="</tr>";
# 	}
# 	if ($table_body != ''){
# 		echo "
# 		<table border=\"1\" class=\"display\" id=\"example2\"><tr><th>Variant</th>";
# 		foreach ($uniprot_header as $field=>$exists){
			
# 			if (preg_match("/accession/",$field)){
# 				echo "<th>Uniprot $field(s)</th><th>Uniprot structure</th><th>Reactome</th>";
# 			}else{
# 				echo "<th>$field</th>";
# 			}
# 		}
# 		echo "</tr>$table_body
# 		</table>
# 		";
# 	}else{
# 		echo "The uniprot file was not found or there were no coding variants in your dataset.";
# 	}
# }else{
# 	echo "The uniprot file was not found or there were no coding variants in your dataset.";
}