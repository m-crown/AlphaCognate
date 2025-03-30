import { useParams } from "react-router-dom";
import { useEffect, useState, useRef } from "react";
import MolStarViewer from "../components/MolStarViewer";
import TransplantTable from "../components/TransplantTable";
import { Grid, Loader, Text } from "@mantine/core";
import FocusButton from "../components/ViewerButton";

type TransplantApiResponse = {
  transplant_data: Array<Transplant>;
  structure_data : Structure;
};

export type Transplant = {
  id: string;
  structure_name: string;
  ligand: string;
  tcs: number;
  struct_asym_id: string; //each transplant is on a separate chain.
};

//this is used in multiple places - abstract away
type Structure = {
  name: string;
  url: string;
  transplanted: boolean;
};

function StructurePage() {
  const { id } = useParams(); // Extract the ID from the URL
  const [transplants, setTransplants] = useState<Transplant[]>([]);
  const [structure, setStructure] = useState<Structure | null>(null);
  const [loading, setLoading] = useState(true);
  //make a url for the structure viewer
  const [structure_url, setStructureUrl] = useState<string | null>(null);
  const viewerInstanceRef = useRef<any>(null);

  useEffect(() => {
    if (id) {
      const fetchTransplantData = async () => {
        const url = new URL(
          '/transplants',
          process.env.NODE_ENV === 'production'
            ? 'http://localhost:8000/' //replace with prod url when it exists
            : 'http://localhost:8000/',
        );
        url.searchParams.append('structure_name', id);
        try {
          const response = await fetch(url.href);
          const json = (await response.json()) as TransplantApiResponse;

          setTransplants(json.transplant_data); // Update state with fetched data
          setStructure(json.structure_data);
          let structure_url = new URL(
            json.structure_data?.url,
            process.env.NODE_ENV === 'production'
              ? 'http://localhost:8080/' //replace with prod url when it exists
              : 'http://localhost:8080/',)
          setStructureUrl(structure_url.href);
          setLoading(false); // Set loading to false after data is fetched
        } catch (error) {
          console.error(error);
          setLoading(false); // Set loading to false even if an error occurs
        }
      };
      fetchTransplantData();
    }
  }, [id]);

  return (
    <Grid>
      <Grid.Col span={6}>
        {loading ? <Loader /> : <TransplantTable data={transplants} viewerInstanceRef={viewerInstanceRef}/>}
        <FocusButton viewerInstanceRef={viewerInstanceRef} />
      </Grid.Col>
      <Grid.Col span={6}>

        <MolStarViewer url={ structure_url || ""} viewerInstanceRef={viewerInstanceRef} /> 
      </Grid.Col>
    </Grid>
  );
}

export default StructurePage;
        //structure?.url