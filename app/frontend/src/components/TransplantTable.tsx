import { useMemo, useState } from 'react';
import {
  MantineReactTable,
  MRT_RowSelectionState,
  useMantineReactTable,
  type MRT_ColumnDef
} from 'mantine-react-table';
import type { TransplantRecord, CognateLigand } from '../pages/StructurePage';

interface TransplantTableProps {
  data: TransplantRecord[]; // Accept data as a prop
  viewerInstanceRef: React.RefObject<any>;
  cognateLigandsData: CognateLigand[]; // New prop for cognate ligands data
  updateCognateLigands: (ligands: CognateLigand[]) => void; // Callback to update cognate ligands in the parent
}

const TransplantTable = ({ data, viewerInstanceRef, cognateLigandsData, updateCognateLigands }: TransplantTableProps) => {
  // table state
  const [currentlyClicked, setCurrentlyClicked] = useState<string | null>(null);
  const [rowSelection, setRowSelection] = useState<MRT_RowSelectionState>({});

  const columns = useMemo<MRT_ColumnDef<TransplantRecord>[]>(
    () => [
      {
        accessorKey: 'transplant.ligand_id',
        header: 'Transplant ID',
      },
      {
        accessorKey: 'transplant.structure_name',
        header: 'Structure Name',
      },
      {
        accessorKey: 'transplant.ligand_id',
        header: 'Ligand',
      },
      {
        accessorKey: 'transplant.tcs',
        header: 'tcs',
      },
      {
        accessorKey: 'transplant.ligand_chain',
        header: 'Chain',
      }
    ],
    [],
  );

  // Filter cognate ligands for the selected ligand
  const filteredCognateLigands = (ligandId: string) => {
    return cognateLigandsData.filter(ligand => ligand.ligand_id === ligandId);
  };

  const table = useMantineReactTable({
    columns,
    data,
    enableMultiRowSelection: false,
    initialState: { showColumnFilters: true },
    mantineTableBodyRowProps: ({ row }) => ({
      onClick: () => {
        const selectedId = row.original.transplant.ligand_chain;
        const alreadySelected = currentlyClicked === selectedId;
        // Update the cognate ligands in the parent component
        if (alreadySelected) {
          viewerInstanceRef.current?.visual.reset({ camera: true });
          setCurrentlyClicked(null);
          setRowSelection({});
        } else {
          viewerInstanceRef.current?.visual.select({
            data: [{ auth_asym_id: selectedId, color: { r: 255, g: 255, b: 0 }, focus: true }],
          });
          setCurrentlyClicked(selectedId);
          setRowSelection({ [row.id]: true });
          console.log('Selected row:', row.original);
          updateCognateLigands(row.original.cognate_ligands);
        }
      },
      selected: !!rowSelection[row.id], // tells MRT the row is "selected"
      sx: (theme) => ({
        cursor: 'pointer',
        backgroundColor: rowSelection[row.id]
          ? theme.colors.yellow[1] // or any selected color you want
          : undefined,
      }),
    }),
    state: {
      rowSelection,
    }
  });

  return (
    <div style={{ maxHeight: '80vh', maxWidth: '100%', overflow: 'auto' }}>
      <MantineReactTable table={table} />
    </div>
  );
};

export default TransplantTable;
