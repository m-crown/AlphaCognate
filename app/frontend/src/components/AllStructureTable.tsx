import { useEffect, useMemo, useState } from 'react';
import {
  MantineReactTable,
  useMantineReactTable,
  type MRT_ColumnDef,
  type MRT_ColumnFiltersState,
  type MRT_PaginationState,
  type MRT_SortingState,
} from 'mantine-react-table';

type StructureApiResponse = {
  data: Array<Structure>;
  meta: {
    totalRowCount: number;
  };
};

type Structure = {
  name: string;
  url: string;
  transplanted: boolean;
  num_transplants: number;
};

//need to have the row action thing here

interface AllStructureTableProps {
  onStructureRowClick: (structure: string) => void; // Callback to pass clicked structure
}


const AllStructureTable: React.FC<AllStructureTableProps> = ({onStructureRowClick}) => {
  //data and fetching state
  const [data, setData] = useState<Structure[]>([]);
  const [isError, setIsError] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const [isRefetching, setIsRefetching] = useState(false);
  const [rowCount, setRowCount] = useState(0);

  //table state
  const [columnFilters, setColumnFilters] = useState<MRT_ColumnFiltersState>(
    [],
  );
  const [globalFilter, setGlobalFilter] = useState('');
  const [sorting, setSorting] = useState<MRT_SortingState>([]);
  const [pagination, setPagination] = useState<MRT_PaginationState>({
    pageIndex: 0,
    pageSize: 10,
  });

  //if you want to avoid useEffect, look at the React Query example instead
  useEffect(() => {
    const fetchData = async () => {
      if (!data.length) {
        setIsLoading(true);
      } else {
        setIsRefetching(true);
      }

      const url = new URL(
        '/structures',
        process.env.NODE_ENV === 'production'
          ? 'http://localhost:8000/' //replace with prod url when it exists
          : 'http://localhost:8000/',
      );
      url.searchParams.set(
        'offset',
        `${pagination.pageIndex * pagination.pageSize}`,
      );
      url.searchParams.set('limit', `${pagination.pageSize}`);
      url.searchParams.set('filters', JSON.stringify(columnFilters ?? []));
      url.searchParams.set('globalFilter', globalFilter ?? '');
      url.searchParams.set('sorting', JSON.stringify(sorting ?? []));

      try {
        const response = await fetch(url.href);
        const json = (await response.json()) as StructureApiResponse;
        setData(json.data);
        setRowCount(json.meta.totalRowCount);
      } catch (error) {
        setIsError(true);
        console.error(error);
        return;
      }
      setIsError(false);
      setIsLoading(false);
      setIsRefetching(false);
    };
    fetchData();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [
    columnFilters, //refetch when column filters change
    globalFilter, //refetch when global filter changes
    pagination.pageIndex, //refetch when page index changes
    pagination.pageSize, //refetch when page size changes
    sorting, //refetch when sorting changes
  ]);

  const columns = useMemo<MRT_ColumnDef<Structure>[]>(
    () => [
      {
        accessorKey: 'name',
        header: 'Structure Name',
      },
      {
        accessorKey: 'url',
        header: 'URL',
      },
      {
        accessorKey: 'transplanted',
        header: 'Transplanted',
      },
      {
        accessorKey: 'num_transplants',
        header: 'Num Transplants',
      }
    ],
    [],
  );

  const table = useMantineReactTable({
    columns,
    data,
    enableRowSelection: true,
    getRowId: (row) => row.id,
    initialState: { showColumnFilters: true },
    manualFiltering: true,
    manualPagination: true,
    manualSorting: true,
    rowCount,
    onColumnFiltersChange: setColumnFilters,
    onGlobalFilterChange: setGlobalFilter,
    onPaginationChange: setPagination,
    onSortingChange: setSorting,
    state: {
      columnFilters,
      globalFilter,
      isLoading,
      pagination,
      showAlertBanner: isError,
      showProgressBars: isRefetching,
      sorting,
    },
    mantineTableBodyRowProps: ({ row }) => {
      return {
          onClick: () => { //this was originally a cell click but changed due to added expansion.
              const selected_structure = row.original.name as string;
              console.log(selected_structure);
              onStructureRowClick(selected_structure); // Pass the gencode ID to the parent
          }
    }},
    mantineToolbarAlertBannerProps: isError
      ? { color: 'red', children: 'Error loading data' }
      : undefined,
  });

  return <div style={{ maxHeight: '80vh', maxWidth: '100%', overflow: 'auto' }}>
      <MantineReactTable table={table} />
    </div>
};

export default AllStructureTable;